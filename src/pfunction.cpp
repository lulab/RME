/*Code for partition function calculation of base pair probabilities.
Copyright 2003,2004,2005,2006 by David H. Mathews

Contributions by Zhi John Lu, 2006.

*/

//NOTE: Two compiler flags exist to modify the behaviour of MB Loops:
//SIMPLEMBLOOP uses a model in which every helix end has a 5' and 3' dangle (Except at sequence ends)
//disablecoax is a flag that disables coaxial stacking, but leaves the rest of the model intact



#include "pfunction.h"
#include "boltzmann.h" //for boltzman
#include <math.h>
#include <cstdlib>

using namespace std;

#undef pfdebugmode  //flag to indicate debugging
//#define pfdebugmode  //flag to indicate debugging

#undef equiout
//#define equiout //flag to indicate equilibria should be written to file

#undef timer

//#define timer //flag to indicate the code execution should be timed

#undef disablecoax
//#define disablecoax



//These have been moved to rna_library
//inline void write(ofstream *out,double *i) {

//	out->write((char *) i,sizeof(*i));
//}

//inline void read(ifstream *out,double *i) {

//	out->read((char *) i,sizeof(*i));
//}

#define maxinter 30   //maximum length of the internal loops
#define maxasym 30  //maximum asymetry in the internal loops


void calculatepfunction(structure* ct,pfdatatable* data, TProgressDialog* update, char* save, bool quickQ, PFPRECISION *Q,
	pfunctionclass *w, pfunctionclass *v, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl,
	pfunctionclass *wcoax, forceclass *fce,PFPRECISION *w5,PFPRECISION *w3,bool *mod, bool *lfce) {


int ip,jp,ii,jj,jpf,jf,bl,ll,dp;
register int i,j,h,d;
int k,p;
register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
int before,after;
register PFPRECISION e;
register int number,maxj;
register PFPRECISION twoscaling,rarray;
PFPRECISION **curE,**prevE,**tempE,**wca;
bool calculatev;


#ifdef timer
	#include <time.h>
	ofstream timeout;
	int seconds;
	char timerstring[100];
	char timelength[10];
	strcpy(timerstring,"time_pf_");
	sprintf(timelength,"%i",ct->numofbases);
	strcat(timerstring,timelength);
	strcat(timerstring,".out");

	timeout.open(timerstring);
	timeout<<time(NULL)<<"\n";
	seconds = time(NULL);
#endif

number = (ct->numofbases);//place the number of bases in a registered integer

curE= new  PFPRECISION *[number+1];
prevE= new PFPRECISION *[number+1];
wca = new PFPRECISION *[number+1];




for (i=0;i<=number;i++) {
	curE[i]= new  PFPRECISION [number+1];
   	prevE[i]= new PFPRECISION [number+1];
	wca[i]=new PFPRECISION [number+1];
	for(j=0;j<=number;j++) {

		wca[i][j] = (PFPRECISION) 0;
		curE[i][j]= (PFPRECISION) 0;
		prevE[i][j]= (PFPRECISION) 0;

	}

}

w5[0] = (PFPRECISION) 1;//initialize the random coil contribution to the partition function
w3[number+1] = (PFPRECISION) 1;

force(ct,fce,lfce);

twoscaling = data->scaling*data->scaling;


//This is the fill routine:

if (quickQ) maxj = number;
else maxj = 2*number-1;

for (h=0;h<=( quickQ?(maxj-1):(maxj-1-minloop) );h++){

	d=(h<=(number-1))?h:(h-number+1);

	if (h==number) {
		for(i=0;i<=number;i++) {
			for(j=0;j<=number;j++) {
				curE[i][j]= (PFPRECISION) 0;
				prevE[i][j]= (PFPRECISION) 0;
			}
		}
	}

	if (((h%10)==0)&&update) update->update((100*h)/(2*ct->numofbases));

	for (i=((h<=(number-1))?1:(2*number-h));i<=((h<=(number-1))?(number-h):number);i++){
		j=i+d;

		//Test to make sure the fragment is large enough to form pairs:
		if (!(j<=(number)&&((j-i)<=minloop))) {



			//First, calculate V(i,j): (The partition function from nucleotides i to j, given that i and j pair)



			//Start the value of V:
			rarray=0.0;


			//Now, test some conditions as to whether V should be evaluated:
			calculatev=true;

			if (ct->templated) {
				//Test whether these nucleotides are allowed to pair because there is a pairing mask in ct ("tem").
				if (i>ct->numofbases) ii = i - ct->numofbases;
				else ii = i;
				if (j>ct->numofbases) jj = j - ct->numofbases;
				else jj = j;
				if (jj<ii) {
					p = jj;
      				jj = ii;
					ii = p;
				}
   				if (!ct->tem[jj][ii]) calculatev = false;
			}


			if (fce->f(i,j)&SINGLE) {
				//i or j is forced single-stranded
				//v->f(i,j) = 0;
				calculatev = false;
			}
			if (fce->f(i,j)&NOPAIR) {
				//i or j is forced into a pair elsewhere
   				//v->f(i,j)= 0;
				calculatev = false;
			}



			if (inc[ct->numseq[i]][ct->numseq[j]]==0) {
				//These are two nucleotides that cannot form a canonical pair
				//v->f(i,j)= 0;
				calculatev = false;
			}


			//force u's into gu pairs, if the user has specified these restraints
			for (ip=0;ip<ct->ngu;ip++) {
				if (ct->gu[ip]==i) {
					if (ct->numseq[j]!=3) {

         				calculatev = false;
					}
				}

				else if (ct->gu[ip]==j) {
       				if (ct->numseq[i]!=3) {

         				calculatev = false;
					}
				}
				else if ((ct->gu[ip]+number)==j) {
       				if (ct->numseq[i]!=3) {

						calculatev = false;
					}
				}

			}



			//now check to make sure that this isn't an isolated pair:
			//	(consider a pair separated by a bulge as not! stacked)

			//before = 0 if a stacked pair cannot form 5' to i
			before =0;
			if ((i>1&&j<(2*number)&&j!=number)) {
				if ((j>number&&((i-j+number)>minloop+2))||j<number) {
					before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
				}
			}


			//after = 0 if a stacked pair cannot form 3' to i
			if ((((j-i)>minloop+2)&&(j<=number)||(j>number+1))&&(i!=number)) {
				after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

			}
			else after = 0;

			//if there are no stackable pairs to i.j then don't allow a pair i,j
			if ((before==0)&&(after==0)) {
				//v->f(i,j)= 0;
				calculatev = false;
			}


			//A large code block for filling V, if all the above conditions pass:
			if (calculatev) {


				//Test to make sure this isn't the end of the sequence
				if (!(i==(number)||j==((number)+1))) {





   					//Perhaps i and j close a hairpin:
					rarray=erg3(i,j,ct,data,fce->f(i,j));

					if ((j-i-1)>=(minloop+2)||j>(number)) {
      					//Perhaps i,j stacks over i+1,j-1
						if (!mod[i]&&!mod[j])  //make sure this is not a site of chemical modification
							rarray+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
						else {
							//allow G-U to be modified or a pair next to a G-U to be modified
							if ((ct->numseq[i]==3&&ct->numseq[j]==4)||(ct->numseq[i]==4&&ct->numseq[j]==3)) {
								rarray+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);

							}
							else if ((ct->numseq[i+1]==3&&ct->numseq[j-1]==4)||(ct->numseq[i+1]==4&&ct->numseq[j-1]==3)) {

								rarray+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);

							}
							else if (i-1>0&&j+1<2*number) {
								if ((ct->numseq[i-1]==3&&ct->numseq[j+1]==4)||(ct->numseq[i-1]==4&&ct->numseq[j+1]==3)) {

									rarray+=erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);

								}

							}

						}
					}



					/* Perhaps i,j closes an interior or bulge loop, search all possibilities
					fill the interior loops' energy rarray first
					calculate the small loop (size<=5) first
					the larger loop is prefilled with curE[dp][i] (after sub2)

					d= j-i, dp= jp-ip (interior loop)
					i<ip<number<jp<j or i<ip<jp<j<number
					*/
					if ((d-1)>=(minloop+3)||j>number)
					for (dp=d-3;dp>=((j>number)?1:minloop+1);dp--) {
						ll=d-dp-2;

						//calculate every ip,jp when ll <=5: 0x1,0x2,0x3,0x4,0x5,1x1,1x2,1x3,1x4,2x2,2x3
						if(ll>=1&&ll<=5) {

							for (ip=i+1;ip<=j-1-dp;ip++){
	  							jp=ip+dp;
								if (inc[ct->numseq[ip]][ct->numseq[jp]])
	  							if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) ) {
   									//using jpf and jf and  instead of jp and j when j,jp>number
									jpf=( (jp<=number)?jp:jp-number);
									jf=((j<=number)?j:j-number);
									rarray+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),
              								fce->f(jpf,jf)) * v->f(ip,jp);
									//i or j is modified
									if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
		  								rarray+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jpf,jf)) *
                  							v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);

      							}

							}
						}
						//when size >=6 and <=30;

						else if (ll>=6&&ll<=maxinter)
						{

							//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is stored in rna_library_inter.cpp
							rarray+=(curE[dp][i]) *erg2ex(i,j,ll,ct,data);

						//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
							for (bl=0;bl<=1;bl++) {
								ip=i+1+bl;
								jp=ip+dp;
								jpf=( (jp<=number)?jp:jp-number);
								jf=((j<=number)?j:j-number);
								if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )

								if (inc[ct->numseq[ip]][ct->numseq[jp]])
								if (abs(ip-i+jp-j)<=maxasym) {
		  							rarray +=(erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jpf,jf))) * (v->f(ip,jp));
									//i or j is modified
		   							if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
										rarray +=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jpf,jf)) *
                  							v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);

								}

								jp=j-1-bl;
								ip=jp-dp;
								jpf=( (jp<=number)?jp:jp-number);
								jf=((j<=number)?j:j-number);
								if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
								if (inc[ct->numseq[ip]][ct->numseq[jp]])
								if (abs(ip-i+jp-j)<=maxasym) {
									rarray += erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jpf,jf)) * v->f(ip,jp);
			   						if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
										rarray+=erg2(i,j,ip,jp,ct,data,fce->f(i,ip),fce->f(jpf,jf)) *
                  							v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);

								}
							}

						}
					}



				}//end of condition to make sure this wasn't the end of the sequence

				//sub1:


				//Perhaps i,j closes a multibranch or exterior loop, enumerate all possibilities

				if (((j-i-1)>=(2*minloop+4))||(j>(number))) {


					//consider the exterior loop closed by i,j
					if (j>number) {
						#ifdef SIMPLEMBLOOP
						//5' and 3' dangling ends
						if (i!=number&&j!=number+1) rarray+= erg4(i,j,i+1,1,ct,data,lfce[i+1])*erg4(i,j,j-1,2,ct,data,lfce[j-1])*
							penalty(i,j,ct,data)*w3[i+1]*w5[j-number-1]*twoscaling;
						//3' dangling end
						else if (i!=number) rarray+= erg4(i,j,i+1,1,ct,data,lfce[i+1])*
							penalty(i,j,ct,data)*w3[i+1]*w5[j-number-1]*twoscaling;
						//5' dangling ends
						else if (j!=number+1) rarray+= erg4(i,j,i+1,1,ct,data,lfce[i+1])*
							penalty(i,j,ct,data)*w3[i+1]*w5[j-number-1]*twoscaling;
						//no dangling ends
						else rarray+= w3[i+1]*w5[j-number-1]*penalty(i,j,ct,data)*twoscaling;

						#else //not simplembloop
         				rarray+= w3[i+1]*w5[j-number-1]*penalty(i,j,ct,data)*twoscaling;


						if (i!=number) rarray+= erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*w3[i+2]*w5[j-number-1]*twoscaling;
						if (j!=(number+1)) rarray+= erg4(i,j,j-1,2,ct,data,lfce[j-1])*penalty(i,j,ct,data)*w3[i+1]*w5[j-number-2]*twoscaling;
						if ((i!=number)&&(j!=(number+1))) {
            				rarray+= data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]*pfchecknp(lfce[i+1],lfce[j-1])*w3[i+2]*
								w5[j-number-2]*penalty(i,j,ct,data)*twoscaling;

						}


						//consider the coaxial stacking of a helix from i to ip onto helix ip+1 or ip+2 to j:
						#ifndef disablecoax //a flag that can turn of coaxial stacking
						//first consider a helix stacking from the 5' sequence fragment:
						for (ip=j-number-minloop-1;ip>0;ip--) {
							//first consider flush stacking
							rarray+=
								w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(j-number-1,ip,ct,data)*
								ergcoaxflushbases(ip,j-number-1,j-number,i,ct,data)*v->f(ip,j-number-1)*twoscaling;

							if ((mod[ip]||mod[j-number-1])) if (j-number-2>0&&notgu(ip,j-number-1,ct)&&!(fce->f(ip,j-number-1)&SINGLE)) {
								if (inc[ct->numseq[ip+1]][ct->numseq[j-number-2]]) {
									rarray+=
										w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(j-number-1,ip,ct,data)*
										ergcoaxflushbases(ip,j-number-1,j-number,i,ct,data)*v->f(ip+1,j-number-2)*
										erg1(ip,j-number-1,ip+1,j-number-2,ct,data)*twoscaling;
								}

							}


							if (j-number-2>0) {
								//now consider an intervening nuc
								if(i<number) {
									rarray+=
										w3[i+2]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip,j-number-2,ct,data)*
										ergcoaxinterbases2(ip,j-number-2,j-number,i,ct,data)*v->f(ip,j-number-2)*twoscaling
										*pfchecknp(lfce[j-number-1],lfce[i+1]);


									if ((mod[ip]||mod[j-number-2])) if (inc[ct->numseq[ip+1]][ct->numseq[j-number-3]]&&notgu(ip,j-number-2,ct)
										&&!(fce->f(ip,j-number-2)&SINGLE)) {
										rarray+=
											w3[i+2]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip,j-number-2,ct,data)*
											ergcoaxinterbases2(ip,j-number-2,j-number,i,ct,data)*v->f(ip+1,j-number-3)*
											erg1(ip,j-number-2,ip+1,j-number-3,ct,data)*twoscaling*pfchecknp(lfce[j-number-1],lfce[i+1]);


									}
								}


								//consider the other possibility for an intervening nuc
								rarray+=
									w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip+1,j-number-2,ct,data)*
									ergcoaxinterbases1(ip+1,j-number-2,j-number,i,ct,data)*v->f(ip+1,j-number-2)*twoscaling
									*pfchecknp(lfce[j-number-1],lfce[ip]);


								if ((mod[ip+1]||mod[j-number-2])) if (inc[ct->numseq[ip+2]][ct->numseq[j-number-3]]&&notgu(ip+1,j-number-2,ct)
									&&!(fce->f(ip+1,j-number-2)&SINGLE)) {
									rarray+=
										w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip+1,j-number-2,ct,data)*
										ergcoaxinterbases1(ip+1,j-number-2,j-number,i,ct,data)*v->f(ip+2,j-number-3)
										*erg1(ip+1,j-number-2,ip+2,j-number-3,ct,data)*twoscaling
										*pfchecknp(lfce[j-number-1],lfce[ip]);
								}


							}


						}

						//now consider a helix stacking from the 3' sequence fragment:
						for (ip=i+minloop+1;ip<=number;ip++) {
							//first consider flush stacking

							rarray+=
								w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip,i+1,ct,data)*
								ergcoaxflushbases(j-number,i,i+1,ip,ct,data)*v->f(i+1,ip)*twoscaling;


							if ((mod[i+1]||mod[ip])) if (inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)
								&&!(fce->f(i+1,ip)&SINGLE)) {

								rarray+=
									w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip,i+1,ct,data)*
									ergcoaxflushbases(j-number,i,i+1,ip,ct,data)*v->f(i+2,ip-1)
									*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling;

							}

							//now consider an intervening nuc
							if (j-number>1) {
								rarray+=
									w3[ip+1]*w5[j-number-2]*penalty(i,j,ct,data)*penalty(ip,i+2,ct,data)*
									ergcoaxinterbases1(j-number,i,i+2,ip,ct,data)*v->f(i+2,ip)*twoscaling
									*pfchecknp(lfce[i+1],lfce[j-number-1]);


								if ((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
									&&!(fce->f(i+2,ip)&SINGLE)) {

									rarray+=
										w3[ip+1]*w5[j-number-2]*penalty(i,j,ct,data)*penalty(ip,i+2,ct,data)*
										ergcoaxinterbases1(j-number,i,i+2,ip,ct,data)*v->f(i+3,ip-1)
										*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling
										*pfchecknp(lfce[i+1],lfce[j-number-1]);

								}
							}


							//consider the other possibility for an intervening nuc
							rarray+=
								w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip-1,i+2,ct,data)*
								ergcoaxinterbases2(j-number,i,i+2,ip-1,ct,data)*v->f(i+2,ip-1)*twoscaling
								*pfchecknp(lfce[i+1],lfce[ip]);

							if ((mod[i+2]||mod[ip-1])) if(inc[ct->numseq[i+3]][ct->numseq[ip-2]]&&notgu(i+2,ip-1,ct)
								&&!(fce->f(i+2,ip-1)&SINGLE)) {

								rarray+=
									w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip-1,i+2,ct,data)*
									ergcoaxinterbases2(j-number,i,i+2,ip-1,ct,data)*v->f(i+3,ip-2)
									*erg1(i+2,ip-1,i+3,ip-2,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip]);


							}



						}
						#endif //ifndef disablecoax
						#endif //simplembloop

					}






					//consider the multiloop closed by i,j
					if ((j-i)>(2*minloop+4)&&i!=number) {
						#ifdef SIMPLEMBLOOP
						if (j-1!=number&&i!=number) {
							//5' and 3' dangling ends
							rarray+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*
								erg4(i,j,j-1,2,ct,data,lfce[j-1])*penalty(i,j,ct,data)*
            					wmb->f(i+1,j-1)* data->eparam[5] * data->eparam[10]*twoscaling;
						}
						else if (j-1!=number) {//i==number
							//5' dangling end
							rarray+=
								erg4(i,j,j-1,2,ct,data,lfce[j-1])*penalty(i,j,ct,data)*
            					wmb->f(i+1,j-1)* data->eparam[5] * data->eparam[10]*twoscaling;
						}
						if (i!=number) {
							//3' dangling ends
							rarray+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*
								penalty(i,j,ct,data)*
            					wmb->f(i+1,j-1)* data->eparam[5] * data->eparam[10]*twoscaling;
						}

						else {
							//no dangling end
							rarray+=wmb->f(i+1,j-1)*
								data->eparam[5]*data->eparam[10]
            					*penalty(i,j,ct,data)*twoscaling;
						}
						#else //!SIMPLEMBLOOP


          				//no dangling ends on i-j pair:
						if (j-1!=number) {
							rarray+=wmb->f(i+1,j-1)*data->eparam[5]*data->eparam[10]
            					*penalty(i,j,ct,data)*twoscaling;


							//i+1 dangles on i-j pair:

							if (i+1!=number) rarray+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*
            					wmb->f(i+2,j-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
						}
						if (j-2!=number) {
							//j-1 dangles
							if (j!=(number+1))rarray+=erg4(i,j,j-1,2,ct,data,lfce[j-1]) * penalty(i,j,ct,data) *
            					wmb->f(i+1,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
							//both i+1 and j-1 dangle
							if ((i+1!=number)&&(j!=(number+1))) {
            					rarray+=data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]*
											pfchecknp(lfce[i+1],lfce[j-1])*
											wmb->f(i+2,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[6]* data->eparam[10]
											*penalty(i,j,ct,data)*twoscaling;
							}
						}





						//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
						#ifndef disablecoax //a flag to turn off coaxial stacking
						for (ip=i+1;(ip<j);ip++) {
							//first consider flush stacking




							//conditions guarantee that the coaxial stacking isn't considering an exterior loop
							//if ((i!=number)/*&&(i+1!=number)*//*&&((j>number)||(ip!=number)&&(ip+1!=number))&&(j-1!=number)*/) {
							if (i!=number&&ip!=number&&j-1!=number) {
								rarray+=penalty(i,j,ct,data)*v->f(i+1,ip)*
									penalty(i+1,ip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(w->f(ip+1,j-1)+wmb->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
									*twoscaling;


								if((mod[i+1]||mod[ip])) if (inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {

									rarray+=penalty(i,j,ct,data)*v->f(i+2,ip-1)*
										penalty(i+1,ip,ct,data)*data->eparam[5]
										*data->eparam[10]*data->eparam[10]*(w->f(ip+1,j-1)+wmb->f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
										*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling;

								}



								//if ((ip<j-1)&&(i+2!=number)) {
								if (ip+2<j-1&&i+1!=number&&ip+1!=number) {
								//now consider an intervening nuc
									if ((ip+2<j-1)/*&&(j>number||ip+2!=number)*/)
									rarray+=penalty(i,j,ct,data)*v->f(i+2,ip)*
										penalty(i+2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w->f(ip+2,j-1)+wmb->f(ip+2,j-1))
										*ergcoaxinterbases2(j,i,i+2,ip,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);

									if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
										&&!(fce->f(i+2,ip)&SINGLE)) {

										rarray+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
											penalty(i+2,ip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(ip+2,j-1)+wmb->f(ip+2,j-1))
											*ergcoaxinterbases2(j,i,i+2,ip,ct,data)
											*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);

									}



									if (ip+1<j-2&&j-2!=number)
									rarray+=penalty(i,j,ct,data)*v->f(i+2,ip)*
										penalty(i+2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w->f(ip+1,j-2)+wmb->f(ip+1,j-2))
										*ergcoaxinterbases1(j,i,i+2,ip,ct,data)*twoscaling
										*pfchecknp(lfce[i+1],lfce[j-1]);

									if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
										&&!(fce->f(i+2,ip)&SINGLE)) {

										rarray+=penalty(i,j,ct,data)*v->f(i+3,ip-1)*
											penalty(i+2,ip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(ip+1,j-2)+wmb->f(ip+1,j-2))
											*ergcoaxinterbases1(j,i,i+2,ip,ct,data)
											*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[j-1]);

									}




								}





							}


						}



						//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
						for (ip=j-1;ip>i;ip--) {


							//conditions guarantee that the coaxial stacking isn't considering an exterior loop
							//if ((i!=number)&&(i+1!=number)&&((j>number)||(ip!=number)&&(ip-1!=number))&&(j-1!=number)) {
							if (j-1!=number&&ip-1!=number&&i!=number) {
								//first consider flush stacking
								rarray+=penalty(i,j,ct,data)*v->f(ip,j-1)*
									penalty(j-1,ip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(w->f(i+1,ip-1)+wmb->f(i+1,ip-1))*ergcoaxflushbases(ip,j-1,j,i,ct,data)
									*twoscaling;


								if((mod[ip]||mod[j-1])) if(inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&notgu(ip,j-1,ct)&&!(fce->f(ip,j-1)&SINGLE)) {
									rarray+=penalty(i,j,ct,data)*v->f(ip+1,j-2)*
										penalty(j-1,ip,ct,data)*data->eparam[5]
										*data->eparam[10]*data->eparam[10]*(w->f(i+1,ip-1)+wmb->f(i+1,ip-1))
										*ergcoaxflushbases(ip,j-1,j,i,ct,data)
										*erg1(ip,j-1,ip+1,j-2,ct,data)*twoscaling;

								}





								if (j-2!=number) {
									//now consider an intervening nuc
									//if ((ip>i+1)&&(j>number||ip-2!=number))
									if (ip-2>i+1&&ip-2!=number) {
										rarray+=penalty(i,j,ct,data)*v->f(ip,j-2)*
											penalty(j-2,ip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(i+1,ip-2)+wmb->f(i+1,ip-2))
											*ergcoaxinterbases1(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);



										if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)&&!(fce->f(ip,j-2)&SINGLE)) {
											rarray+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
												penalty(j-2,ip,ct,data)*data->eparam[5]
												*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
												(w->f(i+1,ip-2)+wmb->f(i+1,ip-2))
												*ergcoaxinterbases1(ip,j-2,j,i,ct,data)
												*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);

										}
									}



									if ((ip-1>i+2)&&i+1!=number) {
										rarray+=penalty(i,j,ct,data)*v->f(ip,j-2)*
											penalty(j-2,ip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(i+2,ip-1)+wmb->f(i+2,ip-1))
											*ergcoaxinterbases2(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);

										if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)
											&&!(fce->f(ip,j-2)&SINGLE)) {
											rarray+=penalty(i,j,ct,data)*v->f(ip+1,j-3)*
												penalty(j-2,ip,ct,data)*data->eparam[5]
												*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
												(w->f(i+2,ip-1)+wmb->f(i+2,ip-1))
												*ergcoaxinterbases2(ip,j-2,j,i,ct,data)
												*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);

										}
									}


								}



							}




						}
						#endif //ifndef disablecoax
						#endif //SIMPLEMBLOOP


						/*if (ct->intermolecular) {

            				//intermolecular, so consider wmb2,
							//don't add the multiloop penalties because this is a exterior loop

            				e[1] = min(e[1],wmb2->f(i+1,j-1) + penalty(i,j,ct,data)+infinity);


            				//i+1 dangles on i-j pair:
            				if (i!=number) e[2] = min(e[2],erg4(i,j,i+1,1,ct,data,lfce[i+1]) + penalty(i,j,ct,data) +
            					wmb2->f(i+2,j-1)+infinity);
            				//j-1 dangles
            				if (j!=(number+1)) e[3] = min(e[3],erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
            					wmb2->f(i+1,j-2)+infinity);
            				//both i+1 and j-1 dangle
            				if ((i!=number)&&(j!=(number+1))) {
            					e[4] = min(e[4],
            					data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
									pfchecknp(lfce[i+1],lfce[j-1]) +
               						wmb2->f(i+2,j-2) + penalty(i,j,ct,data)+infinity);

							}




						}*/


					}



				}

			}//end of the large block for calculating V if all conditions are met
			//sub2:
			//Assign V:
			v->f(i,j) = rarray;

			//Apply constant (an euilibrium constant for formation of i-j pair), if being used:
			if (ct->constant!=NULL) {

				//use ii and jj as the indices for accessing constant array:
				if (i>ct->numofbases) ii = i - ct->numofbases;
				else ii = i;
				if (j>ct->numofbases) jj = j - ct->numofbases;
				else jj = j;
				if (jj<ii) {
					p = jj;
      				jj = ii;
					ii = p;
				}

				v->f(i,j) = v->f(i,j)*ct->constant[jj][ii];

			}



			/*prefill curE[i] and prev[i] for the first two diognals as ll =4 and 5
				As d =10, only fill curE (ll=4, d=10)
				As d =11, only fill prevE (ll=4||5, d=11)
				As d>11, fill curE(ll=4||5, d>11)
				exchange curE and prevE after d >11  as curE[h][i][dp]=curE=curE[h-2][i+1][dp]
				(curE[i][j][ll]=curE[i+1][j-1][ll-2])
			*/


			if ((d-1)>=(minloop+3)||j>number)
			for (dp=d-3;dp>=((j>number)?1:minloop+1);dp--) {
				ll=d-dp-2;
				//calculate every ip>ip+1,jp<jp-1 when ll ==5 ||4
				if (ll==4||ll==5) {
					for (ip=i+2;ip<=j-2-dp;ip++) {
						jp=ip+dp;
						if (inc[ct->numseq[ip]][ct->numseq[jp]])
						if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) ) {
							if(d==( (j>number)?7:10 )||d>( (j>number)?8:11 ))
							//fill the first diagonal of d and first two of ll for every larger d
							{
								curE[dp][i]+=erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp)) * v->f(ip,jp);

								//i or j is modified
								if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
									curE[dp][i]+=erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp)) * v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);

							}


							else if ( d==((j>number)?8:11) )  //fill the second diagonal of d
							{
				  				prevE[dp][i]+=erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp)) * v->f(ip,jp);

				  				if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))

			      					prevE[dp][i]+=erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp))* v->f(ip+1,jp-1)* erg1(ip,jp,ip+1,jp-1,ct,data);
							}



						}

					}
				}

				//when size >=6 and <=30;
				//   else if (ll>=6&&ll<=(data->eparam[7]))
				else if (ll>=6&&ll<=maxinter)
				{
					//calculate minimum curE[dp][i] of 1 x (n-1) for next step's 2 x (n-2)
					ip=i+2;
					jp=ip+dp;
					if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
					if (abs(ip-i+jp-j)<=maxasym)
					if (inc[ct->numseq[ip]][ct->numseq[jp]])
					{
						curE[dp][i] += erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp))*v->f(ip,jp) ;

						if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
							curE[dp][i]+=erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp)) *
			     				v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);


					}

					jp=j-2;
					ip=jp-dp;
					if ( (ip<=number&&jp>number) || (ip<=number&&j<=number) )
					if (abs(ip-i+jp-j)<=maxasym)
					if (inc[ct->numseq[ip]][ct->numseq[jp]])
					{

						curE[dp][i] += erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp))*v->f(ip,jp) ;

						if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
							curE[dp][i] +=erg2in(i,j,ip,jp,ct,data,fce->f(i,ip),
              						fce->f(j,jp))*
			     				v->f(ip+1,jp-1)*erg1(ip,jp,ip+1,jp-1,ct,data);


					}


				}
			}

			//also block propagation of interior loops that contain nucleotides that need to be double-stranded:
			if (lfce[i]||lfce[j]) for (dp=1;dp<=d;dp++) curE[dp][i] = 0.0;


			if (fce->f(i,j)&PAIR)  {//force a pair between i and j
	  			w->f(i,j) = v->f(i,j)*data->eparam[10]*penalty(i,j,ct,data);
				wl->f(i,j) = v->f(i,j)*data->eparam[10]*penalty(i,j,ct,data);
				wmb->f(i,j) = 0.0;
				wmbl->f(i,j) = 0.0;
				wcoax->f(i,j) = 0.0;
				if (j<=number) wca[i][j] = (PFPRECISION) 0;
	  			//goto sub3;
			}
			else {//not a forced pair

				////fill wmb:
				rarray = (PFPRECISION) 0;
				if (((j-i)>(2*minloop+2))||j>number) {


					#ifdef pfdebugmode
						ofstream dump;
						if (i==10&&j==22) {

							dump.open("dump.out");

						}
					#endif



					//also consider the coaxial stacking of two helixes in wv
					e = (PFPRECISION) 0;
					#ifndef SIMPLEMBLOOP
					#ifndef disablecoax //a flag to diable coaxial stacking
					for (ip=i+minloop+1;ip<j-minloop-1;ip++) {
						//first consider flush stacking

						if (ip!=number) {
							rarray+=v->f(i,ip)*v->f(ip+1,j)*penalty(i,ip,ct,data)
								*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data);

							#ifdef pfdebugmode

							if (i==10&&j==22) {

								dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" 11.d.v.1.d.ii.1 rarray+= "<<v->f(i,ip)*v->f(ip+1,j)*penalty(i,ip,ct,data)
								*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)<<" rarray = "<<rarray<<"\n";


							}
							#endif

							if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

								if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
										&&!(fce->f(ip+1,j)&SINGLE)&&!(fce->f(i,ip)&SINGLE)) {

									rarray+=v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+1,j,ip+2,j-1,ct,data);


								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {

									rarray+=v->f(i+1,ip-1)*v->f(ip+1,j)*penalty(i,ip,ct,data)
										*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data);


								}

								if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce->f(ip+1,j)&SINGLE)) {


									rarray+=v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
										*erg1(ip+1,j,ip+2,j-1,ct,data);

								}


							}



							if (ip+1!=number&&j!=number+1) {
								if (!lfce[ip+1]&&!lfce[j]) {
									//now consider an intervening mismatch
									e+=v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data);

							#ifdef pfdebugmode

							if (i==10&&j==22) {

								dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" 11.d.v.1.d.ii.2 earray+= "<<v->f(i,ip)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										<<" earray = "<<e<<"\n";


							}
							#endif

									if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
										if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
											&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
												&&!(fce->f(i,ip)&SINGLE)&&!(fce->f(ip+2,j-1)&SINGLE)) {

											 e+=v->f(i+1,ip-1)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
												*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
												*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+2,j-1,ip+3,j-2,ct,data);


										}

										if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce->f(i,ip)&SINGLE)) {

											e+=v->f(i+1,ip-1)*v->f(ip+2,j-1)*penalty(i,ip,ct,data)
												*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
												*erg1(i,ip,i+1,ip-1,ct,data);


										}

										if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce->f(ip+2,j-1)&SINGLE)) {


											e+=v->f(i,ip)*v->f(ip+3,j-2)*penalty(i,ip,ct,data)
												*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
												*erg1(ip+2,j-1,ip+3,j-2,ct,data);

										}
									}
								}

								if(!lfce[i]&&!lfce[ip+1]&&i!=number) {
									e+=v->f(i+1,ip)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data);

							#ifdef pfdebugmode

							if (i==10&&j==22) {

								dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" 11.d.v.1.d.ii.3 earray+= "<<v->f(i+1,ip)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										<<" earray = "<<e<<"\n";


							}
							#endif
									if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
										if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
											&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
											&&!(fce->f(i+1,ip)&SINGLE)&&!(fce->f(ip+2,j)&SINGLE)	) {

											e+=v->f(i+2,ip-1)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
												*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
												*erg1(i+1,ip,i+2,ip-1,ct,data)*erg1(ip+2,j,ip+3,j-1,ct,data);



										}
										if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce->f(i+1,ip)&SINGLE)) {

											e+=v->f(i+2,ip-1)*v->f(ip+2,j)*penalty(i+1,ip,ct,data)
												*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
												*erg1(i+1,ip,i+2,ip-1,ct,data);


										}

										if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce->f(ip+2,j)&SINGLE)) {


											e+=v->f(i+1,ip)*v->f(ip+3,j-1)*penalty(i+1,ip,ct,data)
												*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
												*erg1(ip+2,j,ip+3,j-1,ct,data);

										}
									}
								}
							}
						}





					}
					#endif //ifndef disablecoax
					#endif //ifndef SIMPLEMBLOOP
					#ifdef pfdebugmode

						if (i==10&&j==22) {

							dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" right before 11.d.v.1.d.iii execution rarray= "<<rarray<<" earray = "<<e<<"\n";
							dump.close();

						}
					#endif

					if (j<=number) wca[i][j] = rarray+e;

					rarray =(rarray+e*data->eparam[6]*data->eparam[6])*data->eparam[10]*data->eparam[10];
					wcoax->f(i,j) = rarray;

					//search for an open bifurcation:
					for (k=i+1;k<j;k++) {
						//e = 0;
						if (k!=number) {
							if (!lfce[i]&&i!=number)
								rarray+=(wl->f(i,k)-wl->f(i+1,k)*data->eparam[6]*data->scaling+wcoax->f(i,k))*(wl->f(k+1,j)+wmbl->f(k+1,j));

							else rarray+=(wl->f(i,k)+wcoax->f(i,k))*(wl->f(k+1,j)+wmbl->f(k+1,j));

						}
         			}




					if (i!=number)
						if (!lfce[i]) rarray+=wmbl->f(i+1,j)*data->eparam[6]*data->scaling;

					wmbl->f(i,j) = rarray;

					wmb->f(i,j) = rarray;
					if (j!=number+1)
						if (!lfce[j]) wmb->f(i,j)+=wmb->f(i,j-1)*data->eparam[6]*data->scaling;


					/*if (ct->intermolecular) {
         				//intermolecular folding:




         			//search for an open bifurcation:
         			for (k=i;k<=j;k++) {

						if (k!=number) wmb2->f(i,j) = min(wmb2->f(i,j),w2->f(i,k)+work2[k+1][jmt]);


         			}



					if (i!=number)
						if (!(fce->f(i,i)&INTER)) wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) );
						else  wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) + data->init - infinity);

					if (j!=number+1)
						if (!(fce->f(j,j)&INTER)) wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1));
						else wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1) +data->init-infinity);



					w2->f(i,j) = min(w2->f(i,j),wmb2->f(i,j) );

				}*/


				}

				//Compute w[i][j]
				if (j>number||(j-i>minloop)) {

					#ifdef SIMPLEMBLOOP

  					//calculate the energy of j stacked onto the pair of i,j-1
					if (j!=N&&i!=1) {
						//5' and 3' dangling ends
     					wl->f(i,j)= v->f(i,j)* data->eparam[10] *
     						erg4(j,i,j+1,1,ct,data,lfce[j+1])*
							erg4(j,i,i-1,2,ct,data,lfce[i-1])*penalty(i,j,ct,data);

						if ((mod[i]||mod[j])) if(inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)&&!(fce->f(i,j)&SINGLE)) {

							wl->f(i,j)+= v->f(i+1,j-1) * data->eparam[10] *
     							erg4(j,i,j+1,1,ct,data,lfce[j+1])*
								erg4(j,i,i-1,2,ct,data,lfce[i-1])*penalty(i,j-1,ct,data)
								*erg1(i,j,i+1,j-1,ct,data);



						}


					}
					if (j!=N) {//i then ==1
						//3' dangling end
     					wl->f(i,j)= v->f(i,j)* data->eparam[10] *
     						erg4(j,i,j+1,1,ct,data,lfce[j+1])*
							penalty(i,j,ct,data);

						if ((mod[i]||mod[j])) if(inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)&&!(fce->f(i,j)&SINGLE)) {

							wl->f(i,j)+= v->f(i+1,j-1) * data->eparam[10] *
     							erg4(j,i,j+1,1,ct,data,lfce[j+1])*
								penalty(i,j-1,ct,data)
								*erg1(i,j,i+1,j-1,ct,data);



						}


					}
					else if (i!=1) {//then j==N
						//5' dangling end
     					wl->f(i,j)= v->f(i,j)* data->eparam[10] *
							erg4(j,i,i-1,2,ct,data,lfce[i-1])*penalty(i,j,ct,data);

						if ((mod[i]||mod[j])) if(inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)&&!(fce->f(i,j)&SINGLE)) {

							wl->f(i,j)+= v->f(i+1,j-1) * data->eparam[10] *
								erg4(j,i,i-1,2,ct,data,lfce[i-1])*penalty(i,j-1,ct,data)
								*erg1(i,j,i+1,j-1,ct,data);



						}


					}

					else {//i==1 and j==N
						wl->f(i,j)= data->eparam[10]*v->f(i,j)*penalty(j,i,ct,data);

						if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

							wl->f(i,j)+= data->eparam[10]*v->f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data);

						}


					}
					#else  //!SIMPLEMBLOOP

					wl->f(i,j)= data->eparam[10]*v->f(i,j)*penalty(j,i,ct,data);

					if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

						wl->f(i,j)+= data->eparam[10]*v->f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data);

					}



					if (i!=number) {
      					//calculate the energy of i stacked onto the pair of i+1,j

						wl->f(i,j)+= v->f(i+1,j)*data->eparam[10]*data->eparam[6]*
         					erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data);

						if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce->f(i+1,j)&SINGLE)) {

							wl->f(i,j)+= v->f(i+2,j-1) * data->eparam[10] *data->eparam[6] *
         						erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)
								*erg1(i+1,j,i+2,j-1,ct,data);


						}

					}
					if (j!=((number)+1)) {
      					//calculate the energy of j stacked onto the pair of i,j-1
						if (j!=1) {
         					wl->f(i,j)+= v->f(i,j-1)* data->eparam[10] * data->eparam[6] *
         						erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data);

							if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce->f(i,j-1)&SINGLE)) {

								wl->f(i,j)+= v->f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         							erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
									*erg1(i,j-1,i+1,j-2,ct,data);



							}


						}
					}
					if ((i!=(number))&&(j!=((number)+1))) {
      					//calculate i and j stacked onto the pair of i+1,j-1
						if (j!=1&&!lfce[i]&&!lfce[j]) {
         					wl->f(i,j)+= v->f(i+1,j-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         						data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								*penalty(j-1,i+1,ct,data);



							if ((mod[i+1]||mod[j-1])) if((j-2>0)&&!(fce->f(i+1,j-1)&SINGLE)) {
								if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

									wl->f(i,j)+= v->f(i+2,j-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         								data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
										*penalty(j-1,i+1,ct,data)*erg1(i+1,j-1,i+2,j-2,ct,data);

								}
							}
						}
					}

					#endif  //SIMPLEMBLOOP


					if (i!=number&&!lfce[i]) {
         				//if (!(fce->f(i,i)&INTER))
					   //add a nuc to an existing loop:
         				wl->f(i,j)+=  wl->f(i+1,j)*data->eparam[6]*data->scaling;
            			//this is for when i represents the center of an intermolecular linker:
						// else e[4] = w->f(i+1,j) + data->eparam[6] + infinity;
					}

					w->f(i,j) = wl->f(i,j);
					if (j!=number+1&&!lfce[j]) {
             			//if (!(fce->f(j,j)&INTER)) {
               			//add a nuc to an existing loop:
               			w->f(i,j)+= w->f(i,j-1) * data->eparam[6]*data->scaling;
					   //}
					   //else e[5] = w->f(i,j-1) + data->eparam[6] + infinity;

					}
				}

				/* if (ct->intermolecular) {

      				//wmb2[i][j%3] = infinity;
      				//keep track of w2:
					for (ii=1;ii<=5;ii++) e[ii] = 2*infinity;


					if (i!=number) {
      					//calculate the energy of i stacked onto the pair of i+1,j

         				e[1] = v->f(i+1,j) +
         					erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data);

						if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]) {

							e[1] = min(e[1],v->f(i+2,j-1) +
         						erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
								+erg1(i+1,j,i+2,j-1,ct,data));

						}



         				e[4] = w2->f(i+1,j);

					}
      				if (j!=((number)+1)) {
      				//calculate the energy of j stacked onto the pair of i,j-1
         				if (j!=1) {
         					e[2] = v->f(i,j-1)   +
         						erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data);

							if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]) {

								e[2] = min(e[2],v->f(i+1,j-2) +
         							erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
									+erg1(i,j-1,i+1,j-2,ct,data));

							}


               				e[5] = w2->f(i,j-1);

         				}
      				}
      				if ((i!=(number))&&(j!=((number)+1))) {
      					//calculate i and j stacked onto the pair of i+1,j-1
         				if (j!=1) {
         					e[3] = v->f(i+1,j-1)   +
         						data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
											[ct->numseq[j]][ct->numseq[i]]
							+pfchecknp(lfce[i+1],lfce[j-1])
               				+penalty(j-1,i+1,ct,data);



							if ((mod[i+1]||mod[j-1])&&inc[ct->numseq[i+2]][ct->numseq[j-2]]) {

								e[3] = min(e[3],v->f(i+2,j-2) +
         							data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
									+pfchecknp(lfce[i+1],lfce[j-1])
									+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data));

							}
         				}
      				}

					e[1] = min(e[1],(v->f(i,j)+penalty(j,i,ct,data)));

					if (mod[i]||mod[j]&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {

						e[1] = min((v->f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data)),e[1]);

					}


      				w2->f(i,j) = min(e[1],e[2]);
      				w2->f(i,j) = min(w2->f(i,j),e[3]);
      				w2->f(i,j) = min(w2->f(i,j),e[4]);
      				w2->f(i,j) = min(w2->f(i,j),e[5]);




				}*/



			}//end of else "not a forced pair"



		}//end condition to mnake sure the fragment is long enough to form pairs

		//sub3:







		//Compute w5[i], the energy of the best folding from 1->i, and
      		//w3[i], the energy of the best folding from i-->numofbases
		if (i==(1)) {


			if (j<=minloop+1) {
				if (lfce[j]) w5[j]= (PFPRECISION) 0;
				else  w5[j] = w5[j-1]*data->scaling;
			}

			else {
      			if (lfce[j]) rarray = (PFPRECISION) 0;

				else rarray = w5[j-1]*data->scaling;



      			for (k=0;k<=(j-4);k++) {

					#ifdef SIMPLEMBLOOP
					//5' and 3' dangling ends
					if (k>0&&j!=N) {
						rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
						*erg4(j,k+1,k,2,ct,data,lfce[k])
						*v->f(k+1,j)*penalty(j,k+1,ct,data);

						if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

							rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
								*erg4(j,k+1,k,2,ct,data,lfce[k])*v->f(k+2,j-1)
								*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
						}
					}
					//5' dangling end
					if (k>0) {//j==N
						rarray+=w5[k]*erg4(j,k+1,k,2,ct,data,lfce[k])
						*v->f(k+1,j)*penalty(j,k+1,ct,data);

						if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

							rarray+=w5[k]*erg4(j,k+1,k,2,ct,data,lfce[k])*v->f(k+2,j-1)
								*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
						}
					}
					//3' dangling end
					if (j!=N) {//k==0
						rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
						*v->f(k+1,j)*penalty(j,k+1,ct,data);

						if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

							rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
								*v->f(k+2,j-1)
								*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
						}
					}
					//no dangling ends
					else {//k==0 and j==N
						rarray+=w5[k]*v->f(k+1,j)*penalty(j,k+1,ct,data);

						if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

							rarray+=w5[k]*v->f(k+2,j-1)
								*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
						}
					}
					#else //!SIMPLEMBLOOP

      				rarray+=w5[k]*v->f(k+1,j)*penalty(j,k+1,ct,data);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

						rarray+=w5[k]*v->f(k+2,j-1)*penalty(j,k+1,ct,data)
							*erg1(k+1,j,k+2,j-1,ct,data);
					}



					rarray+=w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v->f(k+2,j)*penalty(j,k+2,ct,data);

					if((mod[k+2]||mod[j])) if(inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+2,j,ct)
						&&!(fce->f(k+2,j)&SINGLE)) {
						rarray+=w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v->f(k+3,j-1)
							*penalty(j,k+2,ct,data)*erg1(k+2,j,k+3,j-1,ct,data);

					}


         			rarray+=w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+1,j-1)*penalty(j-1,k+1,ct,data);

					if ((mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {

						rarray+=w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+2,j-2)
							*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data);
					}



					rarray+=w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]]
									*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+2,j-1)*
									penalty(j-1,k+2,ct,data);

					if ((mod[k+2]||mod[j-1])) if(inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce->f(k+2,j-1)&SINGLE)) {

						rarray+=w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]]
									*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+3,j-2)*
									penalty(j-1,k+2,ct,data)*erg1(k+2,j-1,k+3,j-2,ct,data);

					}








					rarray+=w5[k]*wca[k+1][j];

					#endif  //SIMPLEMBLOOP
				}


				w5[j] = rarray;

				//check to see if w5 is about to go out of bounds:
				if (w5[j]>PFMAX) {
					rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
					twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
				}
				else if (w5[j]<PFMIN&&w5[j]>0) {
					rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
					twoscaling=twoscaling*SCALEUP*SCALEUP;

				}
			}


			if (j==number) {

				//w3[0] = 0;
				//w3[number+1] = 0;
				for (ii=(number);ii>=(number-minloop);ii--) {    //number+1 ... number-minloop
      				if (lfce[ii]) w3[ii] = (PFPRECISION) 0;
					else w3[ii]=w3[ii+1]*data->scaling;
				}
				//w3[i]=0;
   				for (ii=((number)-minloop-1);ii>=1;ii--) {

      				if (lfce[ii]) rarray = (PFPRECISION) 0;

   					else rarray = w3[ii+1]*data->scaling;



					for (k=((number)+1);k>=(ii+4);k--) {
						#ifdef SIMPLEMBLOOP
						if (ii>1&&k!=N+1) {
							//5' and 3' dangling ends
							rarray+=v->f(ii,k-1)*erg4(k-1,ii,k,1,ct,data,lfce[k])
								*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])
								*penalty(k-1,ii,ct,data)*w3[k];

							if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
								rarray+=v->f(ii+1,k-2)*erg4(k-1,ii,k,1,ct,data,lfce[k])
								*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])*
								penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

							}
						}
						else if (k!=N+1) {//i==1
							//3' dangling end
							rarray+=v->f(ii,k-1)*erg4(k-1,ii,k,1,ct,data,lfce[k])
								*penalty(k-1,ii,ct,data)*w3[k];

							if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
								rarray+=v->f(ii+1,k-2)*erg4(k-1,ii,k,1,ct,data,lfce[k])*
								penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

							}
						}
						else if (ii>1) {//k==N+1
							//5' dangling end
							rarray+=v->f(ii,k-1)
								*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])
								*penalty(k-1,ii,ct,data)*w3[k];

							if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
								rarray+=v->f(ii+1,k-2)
								*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])*
								penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

							}
						}
						else {//ii==1&&k==N+1
							//No dangling ends
							rarray+=v->f(ii,k-1)
								*penalty(k-1,ii,ct,data)*w3[k];

							if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
								rarray+=v->f(ii+1,k-2)*
								penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

							}
						}
						#else //!SIMPLEMBLOOP


      					rarray+=v->f(ii,k-1)*w3[k]*penalty(k-1,ii,ct,data);

						if((mod[ii]||mod[k-1])) if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
							rarray+=v->f(ii+1,k-2)*w3[k]*penalty(k-1,ii,ct,data)*erg1(ii,k-1,ii+1,k-2,ct,data);

						}


						rarray+=v->f(ii+1,k-1)*erg4(k-1,ii+1,ii,2,ct,data,lfce[ii])*penalty(k-1,ii+1,ct,data) * w3[k];

						if((mod[ii+1]||mod[k-1])) if(inc[ct->numseq[ii+2]][ct->numseq[k-2]]&&notgu(ii+1,k-1,ct)&&!(fce->f(ii+1,k-1)&SINGLE)) {

							rarray+=v->f(ii+2,k-2)*erg4(k-1,ii+1,ii,2,ct,data,lfce[ii])*
								penalty(k-1,ii+1,ct,data) *w3[k]*erg1(ii+1,k-1,ii+2,k-2,ct,data);

						}


						rarray+=v->f(ii,k-2)*erg4(k-2,ii,k-1,1,ct,data,lfce[k-1])*penalty(k-2,ii,ct,data)*w3[k];

						if((mod[ii]||mod[k-2]))if(inc[ct->numseq[ii+1]][ct->numseq[k-3]]&&notgu(ii,k-2,ct)&&!(fce->f(ii,k-2)&SINGLE)) {
							rarray+=v->f(ii+1,k-3)*erg4(k-2,ii,k-1,1,ct,data,lfce[k-1])*
								penalty(k-2,ii,ct,data)*w3[k]*erg1(ii,k-2,ii+1,k-3,ct,data);

						}

						if (!lfce[ii]&&!lfce[k-1]) {
							rarray+=v->f(ii+1,k-2)*data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]]
								[ct->numseq[k-1]][ct->numseq[ii]]
								*w3[k]*
								penalty(k-2,ii+1,ct,data);



							if((mod[ii+1]||mod[k-2]))if(inc[ct->numseq[ii+2]][ct->numseq[k-3]]&&notgu(ii+1,k-2,ct)&&!(fce->f(ii+1,k-2)&SINGLE)) {
								rarray+=v->f(ii+2,k-3)*data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]]
									[ct->numseq[k-1]][ct->numseq[ii]]
									*pfchecknp(lfce[k-1],lfce[ii])*w3[k]*
									penalty(k-2,ii+1,ct,data)*erg1(ii+1,k-2,ii+2,k-3,ct,data);


							}
						}

						//also consider coaxial stacking:
						#ifndef disablecoax //a flag to disable coaxial stacking
						for (ip=k+minloop+1;ip<=number+1;ip++) {


							//first consider flush stacking:
							rarray+=v->f(ii,k-1)*v->f(k,ip-1)*w3[ip]*
								penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
								ergcoaxflushbases(ii,k-1,k,ip-1,ct,data);

							if(mod[ii]||mod[k-1]||mod[k]||mod[ip-1]) {

								if ((mod[ii]||mod[k-1])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[ii+1]][ct->numseq[k-2]]
									&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(ii,k-1,ct)&&notgu(k,ip-1,ct)
									&&!(fce->f(ii,k-1)&SINGLE)&&!(fce->f(k,ip-1)&SINGLE)) {

									rarray+=v->f(ii+1,k-2)*v->f(k+1,ip-2)*w3[ip]*
										penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
										*erg1(ii,k-1,ii+1,k-2,ct,data)*erg1(k,ip-1,k+1,ip-2,ct,data);
								}
								if((mod[ii]||mod[k-1])&&inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {

									rarray+=v->f(ii+1,k-2)*v->f(k,ip-1)*w3[ip]*
										penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
										*erg1(ii,k-1,ii+1,k-2,ct,data);

								}

								if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(k,ip-1,ct)&&!(fce->f(k,ip-1)&SINGLE)) {

									rarray+=v->f(ii,k-1)*v->f(k+1,ip-2)*w3[ip]*
										penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
										*erg1(k,ip-1,k+1,ip-2,ct,data);
								}

							}


							//now consider an intervening mismatch:
							if (ii>0&&!lfce[ii]&&!lfce[k-1]) {
								rarray+=v->f(ii+1,k-2)*v->f(k,ip-1)*w3[ip]*
									penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
									ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data);

								if(mod[ii+1]||mod[k-2]||mod[k]||mod[ip-1]){

									if((mod[ii+1]||mod[k-2])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[ii+2]][ct->numseq[k-3]]
										&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(ii+1,k-2,ct)&&notgu(k,ip-1,ct)
										&&!(fce->f(k,ip-1)&SINGLE)&&!(fce->f(ii+1,k-2)&SINGLE)){
										rarray+=v->f(ii+2,k-3)*v->f(k+1,ip-2)*w3[ip]*
											penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
											ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
											*erg1(ii+1,k-2,ii+2,k-3,ct,data)*erg1(k,ip-1,k+1,ip-2,ct,data);

									}

									if((mod[ii+1]||mod[k-2])&&inc[ct->numseq[ii+2]][ct->numseq[k-3]]&&notgu(ii+1,k-2,ct)&&!(fce->f(ii+1,k-2)&SINGLE)) {
										rarray+=v->f(ii+2,k-3)*v->f(k,ip-1)*w3[ip]*
										penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
										*erg1(ii+1,k-2,ii+2,k-3,ct,data);

									}
									if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(k,ip-1,ct)&&!(fce->f(k,ip-1)&SINGLE)) {
										rarray+=v->f(ii+1,k-2)*v->f(k+1,ip-2)*w3[ip]*
											penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
											ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
											*erg1(k,ip-1,k+1,ip-2,ct,data);


									}


								}

							}
							if (!lfce[k-1]&&!lfce[ip-1]) {

								rarray+=v->f(ii,k-2)*v->f(k,ip-2)*w3[ip]*
									penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
									ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data);

								if (mod[ii]||mod[k-2]||mod[k]||mod[ip-2]) {

									if ((mod[ii]||mod[k-2])&&(mod[k]||mod[ip-2])&&inc[ct->numseq[ii+1]][ct->numseq[k-3]]
										&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&notgu(ii,k-2,ct)&&notgu(k,ip-2,ct)
										&&!(fce->f(ii,k-2)&SINGLE)&&!(fce->f(k,ip-2)&SINGLE)) {

										rarray+=v->f(ii+1,k-3)*v->f(k+1,ip-3)*w3[ip]*
											penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
											ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
											*erg1(ii,k-2,ii+1,k-3,ct,data)*erg1(k,ip-2,k+1,ip-3,ct,data);
									}

									if ((mod[ii]||mod[k-2])&&inc[ct->numseq[ii+1]][ct->numseq[k-3]]&&notgu(ii,k-2,ct)&&!(fce->f(ii,k-2)&SINGLE)) {

										rarray+=v->f(ii+1,k-3)*v->f(k,ip-2)*w3[ip]*
											penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
											ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
											*erg1(ii,k-2,ii+1,k-3,ct,data);
									}

									if ((mod[k]||mod[ip-2])&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&notgu(k,ip-2,ct)&&!(fce->f(k,ip-2)&SINGLE)) {

										rarray+=v->f(ii,k-2)*v->f(k+1,ip-3)*w3[ip]*
											penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
											ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
											*erg1(k,ip-2,k+1,ip-3,ct,data);
									}

								}
							}

						}
						#endif //ifndef disablecoax
						#endif //SIMPLEMBLOOP


					}

					w3[ii] = rarray;
					//check to see if w5 is about to go out of bounds:
					if (w3[ii]>PFMAX) {
						rescale(1,number,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
						twoscaling=twoscaling*SCALEDOWN*SCALEDOWN;
					}
					else if (w3[ii]<PFMIN&&w3[ii]>0) {
						rescale(1,number,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
						twoscaling = twoscaling*SCALEUP*SCALEUP;
					}
				}
			}

		}
		//check to see if any of the 2-D arrays are about to go out of bounds
		//(not checking wca[][],curE[][],prevE[][] although they need to be rescaled too)
		if (v->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (w->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wl->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wcoax->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wmb->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wmbl->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (v->f(i,j)<PFMIN&&v->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (w->f(i,j)<PFMIN&&w->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wl->f(i,j)<PFMIN&&wl->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,curE,prevE,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wcoax->f(i,j)<PFMIN&&wcoax->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wmb->f(i,j)<PFMIN&&wmb->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wmbl->f(i,j)<PFMIN&&wmbl->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}

		#ifdef pfdebugmode
		if (twoscaling>PFMAX||twoscaling<PFMIN) {
			//look for positions that will underflow


			ofstream *ufout = new ofstream();
			ufout->open("c:/underflow_alert.out",ios::app);
			*ufout<<"i= "<<i<<" j= "<<j<<" twoscaling= "<<twoscaling<<"\n";
			ufout->close();
			delete ufout;

		}
		#endif//pfdebugmode
	}


	if (d>(j>number?8:11))
	{
		tempE=curE;
		curE=prevE;
		prevE=tempE;
	}
	if(d> (j>number?7:10))
	for (i=((h<=(number-2))?1:(2*number-h-1));i<=((h<=(number-2))?(number-h-1):number);i++)
	for (dp=1;dp<=d-1;dp++)
	{
		if (i<number)
		curE[dp][i]=curE[dp][i+1];
	}




}

for(ii=0;ii<=number;ii++) {
	delete[] wca[ii];
	delete[] curE[ii];
	delete[] prevE[ii];

}
delete[] wca;
delete[] curE;
delete[] prevE;



#ifdef timer
timeout << time(NULL)<<"\n";
timeout << time(NULL) - seconds;
timeout.close();
#endif

//////////////////////////
//output V, W, WMB, and W2V:
#if defined (pfdebugmode)
	ofstream foo;
	foo.open("c:/arrays.out");
	foo << "i" << "\t"<<"j"<<"\t"<<"v->f(i,j)"<<"\t"<<"w->f(i,j)"<<"\t"<<"wmb->f(i,j)\twmbl->f(i,j)\twcoax->f(i,j)"<<"\t"<<"wl->f(i,j)"<<"\t"<<"v->f(j,i+number)"<<"\t"<<"w->f(j,i+number)"<<"\t"<<"wmb->f(j,i+number)"<<"\t"<<"wl->f(j,i+number)"<<"\t"<<"wmbl->f(j,i+numer)\twcoax->f(j,i+number)"<<"\n";
	for (j=1;j<=number;j++) {
		for (i=1;i<=j;i++) {

			foo << i << "\t"<<j<<"\t"<<v->f(i,j)<<"\t"<<w->f(i,j)<<"\t"<<wmb->f(i,j)<<"\t"<<wmbl->f(i,j)<<"\t"<<wcoax->f(i,j)<<"\t"<<wl->f(i,j)<<"\t"<<v->f(j,i+number)<<"\t"<<w->f(j,i+number)<<"\t"<<wmb->f(j,i+number)<<"\t"<<wl->f(j,i+number)<<"\t"<<wmbl->f(j,i+number)<<"\t"<<wcoax->f(j,i+number)<<"\n";

		}
	}

	foo <<"\n\n\n";
	foo << "i" << "\t" << "w5[i]" << "\t" << "w3[i]" << "\n";
	for (i=0;i<=number;i++) foo << i << "\t" << w5[i] << "\t" << w3[i] << "\n";

	foo.close();

#endif

}




//This function cacluates a partition function for the sequence in CT
//If quickQ == true, return the partition function value in Q
//	otherwise, save the partial partition functions in the datafile named save
//If updates on progress are unwanted, set update=NULL
void pfunction(structure* ct,pfdatatable* data, TProgressDialog* update, char* save, bool quickQ, PFPRECISION *Q)
{



int i,j;
bool *lfce,*mod;//[maxbases+1][maxbases+1];
PFPRECISION *w5,*w3,**wca;
int number;



#ifdef equiout
	ofstream kout;
	kout.open("k.out");
	kout << "sequence length = "<<ct->numofbases<<"\n";
	for (i=1;i<=ct->numofbases;i++) {
		kout << tobase(ct->numseq[i]);
		if (i%20==0) kout << "\n";
	}
	kout << "\n";
	kout.close();
#endif





//array *v,*w;
//inc is an array that saves time by showing which bases can pair before
//	erg is called

//number is the number of bases being folded
//v[i][j] is the best energy for subsequence i to j when i and j are paired
//	for i<j<n, it is the interior framgment between nucleotides i and j
//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
//w[i][j] is the best energy for subsequence i to j

number = (ct->numofbases);//place the number of bases in an integer

//scaling is a per nucleotide scale factor for which W and V are divided
//This is necessary to keep the doubles from overflowing:

//scaling = 0.2;//this factor assumes about 1 kcal/mol/base



//allocate space for the v and w arrays:
pfunctionclass w(number);
pfunctionclass v(number);
pfunctionclass wmb(number);
pfunctionclass wl(number);
pfunctionclass wmbl(number);
pfunctionclass wcoax(number);
forceclass fce(number);

if (ct->intermolecular) {
	//take advantage of templating to prevent intramolecular base pairs

	ct->allocatetem();
	for (i=1;i<ct->inter[0];i++) {
		for (j=i+1;j<=ct->inter[2];j++) {

			ct->tem[j][i]=false;

		}
	}
	for (i=ct->inter[2]+1;i<ct->numofbases;i++) {
		for (j=i+1;j<=ct->numofbases;j++) {

			ct->tem[j][i]=false;

		}
	}


}

//This code converts the SHAPE array of data to equilibrium constants.  This is
//needed for the partition function.  NOTE, however, there is no going back so a structure
//that has been used for partition functions cannot then be used to predict a structure
//by free energy minimization. This is a compromise for efficiency, but clearly something
//that could cause a problem for programmers.


if (ct->shaped) {
	for (i=1;i<=2*ct->numofbases;i++) ct->SHAPE[i]=boltzman(ct->SHAPE[i], data->temp);


}

//add a second array for intermolecular folding:

/*if (ct->intermolecular) {
	w2 = new pfunctionclass(number);
	wmb2 = new pfunctionclass(number);




}

else {

	wmb2=NULL;
	w2=NULL;

}*/




lfce = new bool [2*number+1];
mod = new bool [2*number+1];

for (i=0;i<=2*number;i++) {
	lfce[i] = false;
	mod[i] = false;
}


for (i=1;i<=ct->nmod;i++) {

	if (ct->mod[i]!=1&&ct->mod[i]!=ct->numofbases) {
		mod[ct->mod[i]]=true;
		mod[ct->mod[i]+ct->numofbases]=true;
	}
}




w5 = new PFPRECISION [number+1];
w3 = new PFPRECISION [number+2];
//wca = new PFPRECISION *[number+1];



//The next section handles the case where base pairs are not
				//not allowed to form between nucs more distant
				//than ct->maxdistance
if (ct->limitdistance) {

	if (!ct->templated) ct->allocatetem();

	for (j=minloop+2;j<=ct->numofbases;j++) {
		for (i=1;i<j;i++) {
			if (j-i>=ct->maxdistance) ct->tem[j][i]=false;
		}
	}



}



calculatepfunction(ct,data,update,save,quickQ,Q,&w,&v,&wmb,&wl,&wmbl,&wcoax,&fce,w5,w3,mod,lfce);




if (save!=0) {
	writepfsave(save,ct,w5,w3,&v,&w,&wmb,&wl,&wmbl,&wcoax,&fce,mod,lfce,data);
}

if (quickQ) *Q = w5[ct->numofbases];


delete[] lfce;
delete[] mod;




delete[] w5;
delete[] w3;


/*if (ct->intermolecular) {
	delete w2;
	delete wmb2;
}*/



return;
}


//writepfsave writes a save file with partition function data.
void writepfsave(char *filename, structure *ct,
			 PFPRECISION *w5, PFPRECISION *w3,
			 pfunctionclass *v, pfunctionclass *w, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl, pfunctionclass *wcoax,
			 forceclass *fce, bool *mod, bool *lfce, pfdatatable *data) {

	int i,j,k,l,m,n,o,p;
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};

	ofstream sav(filename,ios::binary);

	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors

	//start with a version number of the savefile
	short vers=pfsaveversion;
	write(&sav,&vers); //save a version of the save file

	//start with structure information
	write(&sav,&(ct->numofbases));
	write(&sav,&(ct->intermolecular));
	write(&sav,&data->scaling);

	write(&sav,&(ct->npair));
	for (i=0;i<=ct->npair;i++) {
		write(&sav,&(ct->pair[i][0]));
		write(&sav,&(ct->pair[i][1]));
	}
	for (i=0;i<=ct->numofbases;i++) {

		write(&sav,&(ct->hnumber[i]));
		sav.write(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->numofbases;i++) write(&sav,&(ct->numseq[i]));

	write(&sav,&(ct->ndbl));
	for (i=0;i<=ct->ndbl;i++) write(&sav,&(ct->dbl[i]));


	if (ct->intermolecular) {
		for (i=0;i<3;i++) write(&sav,&(ct->inter[i]));

	}

	write(&sav,&(ct->nnopair));
	for (i=0;i<=ct->nnopair;i++) write(&sav,&(ct->nopair[i]));

	write(&sav,&(ct->nmod));
	for (i=0;i<=ct->nmod;i++) write(&sav,&(ct->mod[i]));

	write(&sav,&(ct->ngu));
	for (i=0;i<=ct->ngu;i++) write(&sav,&(ct->gu[i]));

	write(&sav,ct->ctlabel[1]);

	write(&sav,&(ct->templated));
	if (ct->templated) {
		for (i=0;i<=ct->numofbases;i++) {
			for (j=0;j<=i;j++) write(&sav,&(ct->tem[i][j]));

		}

	}

	//write the SHAPE data (for pseudo-free energy constraints)
	write(&sav,&(ct->shaped));
	if (ct->shaped) {
		for (i=0;i<=2*ct->numofbases;i++) write(&sav,&(ct->SHAPE[i]));
		for (i=0;i<=2*ct->numofbases;i++) write(&sav,&(ct->SHAPEss[i]));

	}


	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct->numofbases;i++) {
		write(&sav,&(w3[i]));
		write(&sav,&(w5[i]));
		for (j=0;j<=ct->numofbases;j++) {
			write(&sav,&(v->dg[i][j]));
			write(&sav,&(w->dg[i][j]));
			write(&sav,&(wmb->dg[i][j]));
			write(&sav,&(wmbl->dg[i][j]));
			write(&sav,&(wl->dg[i][j]));
			write(&sav,&(wcoax->dg[i][j]));
			writesinglechar(&sav,&(fce->dg[i][j]));
			//if (ct->intermolecular) {
			//	write(&sav,&(w2->dg[i][j]));
			//	write(&sav,&(wmb2->dg[i][j]));

			//}


		}


	}


	// I should replace this with functionality to convert .pfs files to .txt files -- rhiju.
	// this is inelegant -- quick way to get at base pair probabilities...
	if ( false ) { //turn this off before CVS checkin!
	  string bpp_file( "bpp.txt" );
	  ofstream bpp_out( bpp_file.c_str() );
	  for (i=1;i<=ct->numofbases;i++) {
	    for (j=1;j<=ct->numofbases;j++) {
	      bpp_out << ' ' << calculateprobability(i,j,v,w5,ct,data,lfce,mod,data->scaling,fce);
	    }
	    bpp_out << endl;
	  }
	  bpp_out.close();
	}


	write(&sav,&(w3[ct->numofbases+1]));
	for (i=0;i<=2*ct->numofbases;i++) {
		write(&sav,&(lfce[i]));
		write(&sav,&(mod[i]));

	}





	//now write the thermodynamic data:
	write(&sav,&(data->temp));
	for (i=0;i<5;i++) write(&sav,&(data->poppen[i]));
	write(&sav,&(data->maxpen));
	for (i=0;i<11;i++) write(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		write(&sav,&(data->inter[i]));
		write(&sav,&(data->bulge[i]));
		write(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					write(&sav,&(data->dangle[i][j][k][l]));
				}
				for (l=0;l<6;l++) {
					write(&sav,&(data->stack[i][j][k][l]));
					write(&sav,&(data->tstkh[i][j][k][l]));
					write(&sav,&(data->tstki[i][j][k][l]));
					write(&sav,&(data->coax[i][j][k][l]));
					write(&sav,&(data->tstackcoax[i][j][k][l]));
					write(&sav,&(data->coaxstack[i][j][k][l]));
					write(&sav,&(data->tstack[i][j][k][l]));
					write(&sav,&(data->tstkm[i][j][k][l]));
					write(&sav,&(data->tstki23[i][j][k][l]));
					write(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							write(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) write(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										write(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}


						}
					}
				}
			}
		}
	}
	write(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		write(&sav,&(data->itloop[i]));
		write(&sav,&(data->tloop[i]));

	}
	write(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		write(&sav,&(data->itriloop[i]));
		write(&sav,&(data->triloop[i]));

	}
	write(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		write(&sav,&(data->ihexaloop[i]));
		write(&sav,&(data->hexaloop[i]));

	}
	write(&sav,&(data->auend));
	write(&sav,&(data->gubonus));
	write(&sav,&(data->cint));
	write(&sav,&(data->cslope));
	write(&sav,&(data->c3));
	write(&sav,&(data->efn2a));
	write(&sav,&(data->efn2b));
	write(&sav,&(data->efn2c));
	write(&sav,&(data->init));
	write(&sav,&(data->mlasym));
	write(&sav,&(data->strain));
	write(&sav,&(data->prelog));
	write(&sav,&(data->singlecbulge));
	write(&sav,&(data->maxintloopsize));



	sav.close();


}



////////////////////////////////////////////////////////////////////////
//pfunctionclass encapsulates the large 2-d arrays of w and v, used by the
//	partition function

      //the constructor allocates the space needed by the arrays
pfunctionclass::pfunctionclass(int size) {
	//zero indicates whether the array should be set to zero as opposed
		//to being set to infinity, it is false by default



	infinite = (PFPRECISION) 0;

    Size = size;
    /*register*/ int i,j;
    dg = new PFPRECISION *[size+1];

	for (i=0;i<=(size);i++)  {
   		dg[i] = new PFPRECISION [size+1];
   	}
    for (i=0;i<=size;i++) {
         for (j=0;j<size+1;j++) {

             dg[i][j] = (PFPRECISION) 0;
         }
    }

}

//the destructor deallocates the space used
pfunctionclass::~pfunctionclass() {


	int i;

    for (i=0;i<=Size;i++) {
        delete[] dg[i];
    }
     delete[] dg;
}

      //f is an integer function that references the correct element of the array
//inline PFPRECISION &pfunctionclass::f(int i, int j) {



//   if (i>j) {
 //       return infinite;
//    }
 //  else if (i>Size) return f(i-Size,j-Size);
 //  else return dg[i][j-i];

//}

//When considering mismatch at the end of a helix, consult this function to check
//	whether the nucs are required to pair
PFPRECISION pfchecknp(bool lfce1,bool lfce2) {

	if (lfce1||lfce2) return 0;
	else return 1;
}



pfdatatable::pfdatatable() {
}

pfdatatable::pfdatatable(datatable *data, PFPRECISION Scaling, PFPRECISION  Temp) {
	//the partition function datatable needs to be initialized from the datatable
	short i,j,k,l,m,n,o,p;

	scaling = Scaling;

	//store the temperature in the pfdatatable
	temp = Temp;

	for (i=1;i<5;i++) poppen[i]=boltzman(data->poppen[i],temp);
	maxpen = boltzman(data->maxpen,temp);
	for (i=1;i<11;i++) eparam[i] = boltzman(data->eparam[i],temp);
	maxintloopsize = data->eparam[7];
	for (i=1;i<31;i++) {
		inter[i] = boltzman(data->inter[i],temp)*pow(scaling,i+2);
		bulge[i] = boltzman(data->bulge[i],temp)*pow(scaling,i+2);
		hairpin[i] = boltzman(data->hairpin[i],temp)*pow(scaling,i+2);

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=1;l<3;l++) {
					#ifdef SIMPLEMBLOOP
					//In the case of simple multibranch loops, dangles no longer
						//occupy a position in the sequence and do not need a scaling
						//factor
					dangle[i][j][k][l] = boltzman(data->dangle[i][j][k][l],temp);
					#else //!SIMPLEMBLOOP
					dangle[i][j][k][l] = boltzman(data->dangle[i][j][k][l],temp)*scaling;

					#endif
				}
				for (l=0;l<6;l++) {
					stack[i][j][k][l]=boltzman(data->stack[i][j][k][l],temp)*pow(scaling,2);
					tstkh[i][j][k][l]=boltzman(data->tstkh[i][j][k][l],temp);
					tstki[i][j][k][l]=boltzman(data->tstki[i][j][k][l],temp);
					coax[i][j][k][l]=boltzman(data->coax[i][j][k][l],temp);
					tstackcoax[i][j][k][l]=boltzman(data->tstackcoax[i][j][k][l],temp)*pow(scaling,2);
					coaxstack[i][j][k][l] = boltzman(data->coaxstack[i][j][k][l],temp);
					tstack[i][j][k][l]=boltzman(data->tstack[i][j][k][l],temp)*pow(scaling,2);
					tstkm[i][j][k][l]=boltzman(data->tstkm[i][j][k][l],temp)*pow(scaling,2);
					tstki23[i][j][k][l]=boltzman(data->tstki23[i][j][k][l],temp);
					tstki1n[i][j][k][l]=boltzman(data->tstki1n[i][j][k][l],temp);
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							iloop11[i][j][k][l][m][n]=boltzman(data->iloop11[i][j][k][l][m][n],temp)*pow(scaling,4);
							for (o=0;o<6;o++) {
								iloop21[i][j][k][l][m][n][o]=boltzman(data->iloop21[i][j][k][l][m][n][o],temp)*pow(scaling,5);
								for (p=0;p<6;p++) {
									iloop22[i][j][k][l][m][n][o][p]=boltzman(data->iloop22[i][j][k][l][m][n][o][p],temp)*pow(scaling,6);
								}
							}


						}
					}
				}
			}
		}
	}
	numoftloops = data->numoftloops;
	for (i=1;i<=data->numoftloops;i++) {
		itloop[i]=data->tloop[i][0];
		tloop[i] = boltzman(data->tloop[i][1],temp)*pow(scaling,6);

	}
	numoftriloops=data->numoftriloops;
	for (i=1;i<=data->numoftriloops;i++) {
		itriloop[i] = data->triloop[i][0];
		triloop[i] = boltzman(data->triloop[i][1],temp)*pow(scaling,5);
	}
	numofhexaloops=data->numofhexaloops;
	for (i=1;i<=data->numofhexaloops;i++) {
		ihexaloop[i]=data->hexaloop[i][0];
		hexaloop[i] = boltzman(data->hexaloop[i][1],temp)*pow(scaling,8);

	}
	auend = boltzman(data->auend,temp);
	gubonus = boltzman(data->gubonus,temp);
	cint = boltzman(data->cint,temp);
	cslope = boltzman(data->cslope,temp);
	c3 = boltzman(data->c3,temp);
	efn2a = boltzman(data->efn2a,temp);
	efn2b = boltzman(data->efn2b,temp);
	efn2c = boltzman(data->efn2c,temp);
	init = boltzman(data->init,temp);
	mlasym = boltzman(data->mlasym,temp);
	strain = boltzman(data->strain,temp);
	prelog = data->prelog/conversionfactor;
	singlecbulge = boltzman(data->singlecbulge,temp);





}





//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

PFPRECISION erg1(int i,int j,int ip,int jp,structure *ct, pfdatatable *data)
{

		PFPRECISION energy;

		 if ((i==(ct->numofbases))||(j==((ct->numofbases)+1))) {
      		//this is not allowed because n and n+1 are not cavalently attached
			energy = (PFPRECISION) 0;
		}
		else {
      		energy = data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])]*data->eparam[1];

				if (ct->shaped) {
					energy=energy*ct->SHAPE[i];
					energy=energy*ct->SHAPE[j];
					energy=energy*ct->SHAPE[ip];
					energy=energy*ct->SHAPE[jp];
				}

				if ( ct->experimentalPairBonusExists ) {

				    energy = energy * ct->EX[i][j] * ct->EX[ip][jp];

				}
		}

		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k base pair stack ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}


//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,
	char a, char b)
{

	int size,size1,size2, lopsid, count;
	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   	if (((i<=(ct->numofbases))&&(ip>(ct->numofbases)))||((
      	jp<=(ct->numofbases))&&(j>(ct->numofbases)))) {
         //A loop cannot contain the ends of the sequence

         return 0;
      }



      size1 = ip-i-1;
		size2 = j - jp - 1;

	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return 0;//the loop contains a nuc that
      		//should be double stranded



	}


      //a typical internal or bulge loop:
		//size1 = ip-i-1;
		//size2 = j - jp - 1;
		if (size1==0||size2==0) {//bulge loop


			size = size1+size2;
			if (size==1) {
				count = 1;
				energy = data->stack[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[ip]][ct->numseq[jp]]
						*data->bulge[size]*data->eparam[2]/(data->scaling*data->scaling);
				if (size1==1)  {

					//count the number of alternative bulges that exist:

					//k = i;
					//while (ct->numseq[k]==ct->numseq[i+1]) {
					//	count++;
					//	k--;
					//}
					//k=ip;
					//while(ct->numseq[k]==ct->numseq[i+1]) {
					//	count++;
					//	k++;
					//}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[i+1]==2&&(ct->numseq[i+2]==2||ct->numseq[i]==2)) energy= energy*data->singlecbulge;

				}

				else {
					//size2 == 1

					//count the number of alternative bulges that exist:

					//k = jp;
					//while (ct->numseq[k]==ct->numseq[jp+1]) {
					//	count++;
					//	k--;
					//}
					//k=j;
					//while(ct->numseq[k]==ct->numseq[jp+1]) {
					//	count++;
					//	k++;
					//}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[j-1]==2&&(ct->numseq[j-2]==2||ct->numseq[j]==2)) energy=energy*data->singlecbulge;

				}
				//do not apply a correction for the number of equivalent states because
					//the bulge can move to adjacent sites in the partition function calc
				//energy-= (int) (rt*conversionfactor* log ((double) count));
			}
			else if (size>30) {

				loginc = ((data->prelog)*log(PFPRECISION ((size)/30.0)));
				energy = data->bulge[30]*exp(-loginc/(RKC*data->temp))*data->eparam[2];
				energy = energy*penalty(i,j,ct,data)*penalty(jp,ip,ct,data)*pow(data->scaling,size-30);

			}
			else {
         		energy = data->bulge[size]*data->eparam[2];
				energy = energy*penalty(i,j,ct,data)*penalty(jp,ip,ct,data);
			}
		}
		else {//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {

				loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));
				if (size1==1||size2==1) {
            		energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp)) *data->eparam[3]*
						max(data->maxpen,
						pow(data->poppen[min(2,min(size1,size2))],lopsid))
						*pow(data->scaling,size-30);

				}

				else {
					energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp))*data->eparam[3] *
						max(data->maxpen,
						pow(data->poppen[min(2,min(size1,size2))],lopsid))
						*pow(data->scaling,size-30);
				}
			}
			else if ((size1==2)&&(size2==2))//2x2 internal loop
			    energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[jp]]
					[ct->numseq[i+1]][ct->numseq[i+2]]
					[ct->numseq[j-1]][ct->numseq[j-2]];


			else if ((size1==1)&&(size2==2)) {//2x1 internal loop
				energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
					[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];


			}
			else if ((size1==2)&&(size2==1)) {//1x2 internal loop
				energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
					[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];

			}

			else if (size==2) //a single mismatch

				energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];
			else if (size1==1||size2==1) { //this loop is lopsided
         	//this is a 1xn loop:
				energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
						data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
			}


			else if ((size1==2&&size2==3)||(size1==3&&size2==2)) {
			//this is a 2x3 loop
				energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
					data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));


			}
			else {



         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
			}
		}
		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k internal/bulge ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}

//calculate the energy of the intorior part of a internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2in(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,char a,char b)
{

	int size,size1,size2, lopsid;
	PFPRECISION energy;

	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return 0;//the loop contains a nuc that
      		//should be double stranded



	}

   	if (((i<=(ct->numofbases))&&(ip>(ct->numofbases)))||((
      	jp<=(ct->numofbases))&&(j>(ct->numofbases)))) {
         //A loop cannot contain the ends of the sequence

         return 0;
      }
      size1 = ip-i-1;
	  size2 = j - jp - 1;
      size = size1 + size2;
	  lopsid = abs(size1-size2);
	  if (size1 != 0 && size2 !=0)
	 	{
        energy = 	data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					 data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
		}

		return energy;
}


//calculate the energy of the exterior part of internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2ex(int i,int j,int size, structure *ct, pfdatatable *data)
{

	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
				if (size>30) {

				loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp))*pow(data->scaling,size-30);
			}
						else
         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
						 data->inter[size]	;
		return energy;
}




//calculate the energy of a hairpin loop:
#ifndef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#else
PFPRECISION erg3indirect(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#endif
int size,count,key,k;
PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/

	if ((i<=(ct->numofbases))&&(j>(ct->numofbases))) {
      	//A hairpin cannot contain the ends of the sequence

         return 0;
    }

	if (dbl&DUBLE) return 0;//the loop contains a base that should be
      										//double stranded

    else if (dbl&INTER) {//intermolecular interaction
      	//intermolecular "hairpin" free energy is that of intermolecular
         //	initiation times the stacked mismatch

         energy =  data->tstack[ct->numseq[i]][ct->numseq[j]]
         	[ct->numseq[i+1]][ct->numseq[j-1]]*penalty(i,j,ct,data)*pow(data->scaling,j-i-3);


		//Add other states that can exist, i.e. 5' and 3' dangles:
		if (ct->numseq[i+1]!=5&&ct->numseq[j-1]!=5) {
			//Add a 3' dangling end

			energy+=erg4(i,j,i+1,1,ct,data,false)*pow(data->scaling,j-i-2)*penalty(i,j,ct,data);


			//Add a 5' dangling end

			energy+=erg4(i,j,j-1,2,ct,data,false)*pow(data->scaling,j-i-2)*penalty(i,j,ct,data);


			//Add the case where nothing stacks

			energy+=pow(data->scaling,j-i-1)*penalty(i,j,ct,data);

		}
		else if (ct->numseq[i+1]!=5||ct->numseq[j-1]!=5) {

			//Add the case where nothing stacks

			energy+=pow(data->scaling,j-i-1)*penalty(i,j,ct,data);

		}


         return data->init*energy;
    }




		size = j-i-1;



		if (size>30) {

			loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));



			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[30]*exp(-loginc/(RKC*data->temp))*data->eparam[4]*pow(data->scaling,size-30);
		}
		else if (size<3) {
      		energy = data->hairpin[size]*data->eparam[4];
				if (ct->numseq[i]==4||ct->numseq[j]==4) energy = energy*exp(-.6/(RKC*data->temp));
		}
		else if (size==4) {

			key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftloops;count++) {
					if (key==data->itloop[count]) return data->tloop[count];
			}
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}
		else if (size==3) {

			key = (ct->numseq[j])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftriloops;count++) {
				if (key==data->itriloop[count]) return data->triloop[count];
			}

			energy =	data->hairpin[size] * data->eparam[4]
         	*penalty(i,j,ct,data);
		}
		else if (size==6) {
			key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125
				+ (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numofhexaloops;count++) {
				if (key==data->ihexaloop[count]) return data->hexaloop[count];
			}

			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}

		else {
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}




		//check for GU closeure preceded by GG
      if (ct->numseq[i]==3&&ct->numseq[j]==4) {
      	if ((i>2&&i<ct->numofbases)||(i>ct->numofbases+2))
       		if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {

         		energy = energy * data->gubonus;



         	}
      }

      //check for an oligo-c loop

      for (k=1;(k<=size);k++) {
       	if (ct->numseq[i+k] != 2) return energy;//this is not an oligo-C loop
      }
      //this is a poly c loop so penalize
      if (size==3) return (energy *data->c3);
      else return (energy * data->cint * pow(data->cslope,size));


}

#ifdef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
	PFPRECISION energy;

	energy = erg3indirect(i,j,ct, data,dbl,temp);
	ofstream kout;
	kout.open("k.out",ofstream::app);
	kout << "k hairpin ("<<i<<","<<j<<") ="<<energy<< "\n";
	kout.close();
	return energy;
}
#endif

//calculate the energy of a dangling end:
PFPRECISION erg4(int i,int j,int ip,int jp,structure *ct, pfdatatable *data, bool lfce)
{

//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle



      if (lfce) return 0;//stacked nuc should be double stranded

	  //commented out 11/8/99
      //if (ip==5) return 0;//dangling nuc is an intermolecular linker
	    #ifdef equiout
			ofstream kout;
			kout.open("k.out",ofstream::app);
			kout << "k dangle ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp]<< "\n";
			kout.close();
		#endif

		return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];

}

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty(int i,int j,structure* ct, pfdatatable *data) {

	#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		if (ct->numseq[i]==4||ct->numseq[j]==4) kout << "k penalty ("<<i<<","<<j<<") ="<<data->auend<< "\n";
		else kout << "k penalty ("<<i<<","<<j<<") ="<<1.0<< "\n";
		kout.close();
	#endif

	if (ct->numseq[i]==4||ct->numseq[j]==4)
   	return data->auend;
	else return 1;//no end penalty




}

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty2(int i,int j, pfdatatable *data) {


   if (i==4||j==4)
   	return data->auend;
   else return 1;//no end penalty


}



//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//#ifdef equiout
//PFPRECISION ergcoaxindirect(int i, int j, int ip, int jp, int k, structure *ct, pfdatatable *data) {
//#else
PFPRECISION ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif

		//flush stacking

		return data->coax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]];





}

PFPRECISION ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
		//coaxial stacking with an intervening mismatch
		//(k==i-1)

			return data->tstackcoax[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[i-1]] *
				data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]]
				[ct->numseq[ip]][ct->numseq[jp]];


}

PFPRECISION ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
	//coaxial stacking with an intervening mismatch
	/*(k==jp+1) {*/

			return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]]
				[ct->numseq[jp+1]][ct->numseq[ip-1]] *
				data->coaxstack[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[jp+1]];

}

//#ifdef equiout
//PFPRECISION ergcoax(int i, int j, int ip, int jp, int k, structure *ct, pfdatatable *data) {
//PFPRECISION energy;

//	energy = ergcoaxindirect(i,j,ip,jp,k,ct,data);

//	ofstream kout;

//	if (((i<j&&j<ip&&ip<jp)||(i>j&&j>ip&&ip>jp))&&k!=i&&k!=j&&k!=ip&&k!=jp) {
//		kout.open("k.out",ofstream::app);
//		kout << "k coax ("<<i<<","<<j<<","<<ip<<","<<jp<<","<<k<<") ="<<energy<< "\n";
//		kout.close();
//	}
//	return energy;
//}
//#endif


void readpfsave(const char *filename, structure *ct,
			 PFPRECISION *w5, PFPRECISION *w3,
			 pfunctionclass *v, pfunctionclass *w, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl, pfunctionclass *wcoax,
			 forceclass *fce, PFPRECISION *scaling, bool *mod, bool *lfce, pfdatatable *data) {
	 int i,j,k,l,m,n,o,p;
	ifstream sav(filename,ios::binary);
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};//a mask array indicating the identity of canonical pairs
	short vers;

	//read the save file

	//read the file version first
	read(&sav,&vers);

	//start with structure information
	read(&sav,&(ct->numofbases));
	read(&sav,&(ct->intermolecular));
	read(&sav,scaling);

	data->scaling=*scaling;

	read(&sav,&(ct->npair));
	for (i=0;i<=ct->npair;i++) {
		read(&sav,&(ct->pair[i][0]));
		read(&sav,&(ct->pair[i][1]));
	}
	for (i=0;i<=ct->numofbases;i++) {

		read(&sav,&(ct->hnumber[i]));
		sav.read(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->numofbases;i++) read(&sav,&(ct->numseq[i]));

	read(&sav,&(ct->ndbl));
	for (i=0;i<=ct->ndbl;i++) read(&sav,&(ct->dbl[i]));


	if (ct->intermolecular) {
		for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));

	}

	read(&sav,&(ct->nnopair));
	for (i=0;i<=ct->nnopair;i++) read(&sav,&(ct->nopair[i]));

	read(&sav,&(ct->nmod));
	for (i=0;i<=ct->nmod;i++) read(&sav,&(ct->mod[i]));

	read(&sav,&(ct->ngu));
	for (i=0;i<=ct->ngu;i++) read(&sav,&(ct->gu[i]));

	read(&sav,ct->ctlabel[1]);

	read(&sav,&(ct->templated));
	if (ct->templated) {
		ct->allocatetem();
		for (i=0;i<=ct->numofbases;i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct->tem[i][j]));

		}

	}

	read(&sav,&(ct->shaped));
	if (ct->shaped) {
		ct->SHAPE = new double [2*ct->numofbases+1];
		for (i=0;i<=2*ct->numofbases;i++) read(&sav,&(ct->SHAPE[i]));
		ct->SHAPEss = new double [2 * ct->numofbases + 1];
		for (i=0;i<=2*ct->numofbases;i++) read(&sav,&(ct->SHAPEss[i]));

	}


	//now read the array class data for v, w, and wmb:
	for (i=0;i<=ct->numofbases;i++) {
		read(&sav,&(w3[i]));
		read(&sav,&(w5[i]));
		for (j=0;j<=ct->numofbases;j++) {
			read(&sav,&(v->dg[i][j]));
			read(&sav,&(w->dg[i][j]));
			read(&sav,&(wmb->dg[i][j]));
			read(&sav,&(wmbl->dg[i][j]));
			read(&sav,&(wl->dg[i][j]));
			read(&sav,&(wcoax->dg[i][j]));
			readsinglechar(&sav,&(fce->dg[i][j]));
			//if (ct->intermolecular) {
			//	read(&sav,&(w2->dg[i][j]));
			//	read(&sav,&(wmb2->dg[i][j]));

			//}


		}


	}

	read(&sav,&(w3[ct->numofbases+1]));
	for (i=0;i<=2*ct->numofbases;i++) {
		read(&sav,&(lfce[i]));
		read(&sav,&(mod[i]));

	}





	//now read the thermodynamic data:
	read(&sav,&(data->temp));
	for (i=0;i<5;i++) read(&sav,&(data->poppen[i]));
	read(&sav,&(data->maxpen));
	for (i=0;i<11;i++) read(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		read(&sav,&(data->inter[i]));
		read(&sav,&(data->bulge[i]));
		read(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					read(&sav,&(data->dangle[i][j][k][l]));
				}
				for (l=0;l<6;l++) {
					read(&sav,&(data->stack[i][j][k][l]));
					read(&sav,&(data->tstkh[i][j][k][l]));
					read(&sav,&(data->tstki[i][j][k][l]));
					read(&sav,&(data->coax[i][j][k][l]));
					read(&sav,&(data->tstackcoax[i][j][k][l]));
					read(&sav,&(data->coaxstack[i][j][k][l]));
					read(&sav,&(data->tstack[i][j][k][l]));
					read(&sav,&(data->tstkm[i][j][k][l]));
					read(&sav,&(data->tstki23[i][j][k][l]));
					read(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							read(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) read(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										read(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}


						}
					}
				}
			}
		}
	}
	read(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		read(&sav,&(data->itloop[i]));
		read(&sav,&(data->tloop[i]));

	}
	read(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		read(&sav,&(data->itriloop[i]));
		read(&sav,&(data->triloop[i]));

	}
	read(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		read(&sav,&(data->ihexaloop[i]));
		read(&sav,&(data->hexaloop[i]));

	}
	read(&sav,&(data->auend));
	read(&sav,&(data->gubonus));
	read(&sav,&(data->cint));
	read(&sav,&(data->cslope));
	read(&sav,&(data->c3));
	read(&sav,&(data->efn2a));
	read(&sav,&(data->efn2b));
	read(&sav,&(data->efn2c));
	read(&sav,&(data->init));
	read(&sav,&(data->mlasym));
	read(&sav,&(data->strain));
	read(&sav,&(data->prelog));
	read(&sav,&(data->singlecbulge));
	read(&sav,&(data->maxintloopsize));



	sav.close();
}

//return the pairing probability of the i=j pair, where i<j.
PFPRECISION calculateprobability(int i, int j, pfunctionclass *v, PFPRECISION *w5, structure *ct, pfdatatable *data, bool *lfce, bool *mod, PFPRECISION scaling, forceclass *fce) {
	PFPRECISION interior, exterior;
	short before,after;
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	bool adjacentgu;

	if (!mod[i]&&!mod[j]) {
		if (ct->constant==NULL) return (v->f(i,j)*v->f(j,i+ct->numofbases))/(w5[ct->numofbases]*scaling*scaling);
		else {
			//constant is being used.
			if (ct->constant[j][i]<EPSILON) return 0.0;
			return (v->f(i,j)*v->f(j,i+ct->numofbases))/(w5[ct->numofbases]*scaling*scaling*ct->constant[j][i]);
		}
	}
	else {
		if (!(fce->f(i,j)&SINGLE)) {
			before =0;
			if ((i>1&&j<(2*ct->numofbases)&&j!=ct->numofbases)) {
				if ((j>ct->numofbases&&((i-j+ct->numofbases)>minloop+2))||j<ct->numofbases) {
					before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
				}
			}


			//after = 0 if a stacked pair cannot form 3' to i
			if ((((j-i)>minloop+2)&&(j<=ct->numofbases)||(j>ct->numofbases+1))&&(i!=ct->numofbases)) {
				after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

			}
			else after = 0;

			adjacentgu = false;
			//check for preceding or following GU or whether the pair itself is gu
			if (ct->numseq[i+1]==3&&ct->numseq[j-1]==4) adjacentgu = true;
			else if (ct->numseq[i+1]==4&&ct->numseq[j-1]==3) adjacentgu = true;
			else if (ct->numseq[i]==3&&ct->numseq[j]==4) adjacentgu = true;
			else if (ct->numseq[i]==4&&ct->numseq[j]==3) adjacentgu = true;
			else if (i-1>0&&j+1<=ct->numofbases) {
				if (ct->numseq[i-1]==3&&ct->numseq[j+1]==4) adjacentgu = true;
				else if (ct->numseq[i-1]==4&&ct->numseq[j+1]==3) adjacentgu = true;

			}

			//if there are no stackable pairs to i.j then don't allow a pair i,j
			if ((before!=0)||(after!=0)) {

				if (i+1<j-1&&!adjacentgu) interior = erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
				else interior = (PFPRECISION) 0;
				if (j+1<=ct->numofbases&&!adjacentgu) exterior = erg1(j,i+ct->numofbases,j+1,i+ct->numofbases-1,ct,data)*v->f(j+1,i+ct->numofbases-1);
				else exterior = (PFPRECISION) 0;
				return ((v->f(i,j)+interior)*(v->f(j,i+ct->numofbases)+exterior)-interior*exterior)/(w5[ct->numofbases]*scaling*scaling);
			}
			else return 0;


		}
		else return 0;

	}

}

//function to rescale all arrays when partition function calculation is headed out of bounds
void rescale(int i, int j, structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
			 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION **wca, PFPRECISION **curE, PFPRECISION **prevE, PFPRECISION rescalefactor) {
	int h,d,dp,ii,jj,hh,nucs,index,lowi,highi,number;
	double multiplier;
	number=ct->numofbases;
	d=j-i;
	h=(j<=number)?d:(d+number-1);

#ifdef pfdebugmode
	ofstream ufout;
	ufout.open("c:/underflow_alert.out",ios::app);
	ufout<<"rescale factor = "<<rescalefactor<<"\n";
	ufout.close();
#endif
	//rescale v,w,wl,wcoax,wmb,wmbl,wca
	for (hh=0;hh<=h;hh++){
		//determine the boundary of ii
		if (h < number) {
			if (hh == h) {
				highi=i;
				lowi=1;
			}
			else {
				highi=number-hh;
				lowi=1;
			}
		}
		else {
			if (hh==h) {
				highi=i;
				lowi= 2*number-hh;
			}
			else if(hh >= number) {
				highi=number;
				lowi=2*number-hh;
			}
			else if(hh<number) {
				highi=number-hh;
				lowi=1;
			}
		}
		for (ii=lowi;ii<=highi;ii++){
			jj=(hh<=(number-1))?hh+ii:(hh+ii-number+1);
			nucs = jj-ii+1; //this work even jj> number
			multiplier = pow(rescalefactor,(PFPRECISION) nucs);


#ifdef pfdebugmode
			//look for positions that will underflow

			if (multiplier<0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}



			}

			if (multiplier>0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}



			}

#endif //debug
			//rescale v,w,wl,wcoax,wmb,wmbl,wca
			v->f(ii,jj)=v->f(ii,jj)*multiplier;
			w->f(ii,jj)=w->f(ii,jj)*multiplier;
			wl->f(ii,jj)=wl->f(ii,jj)*multiplier;
			wcoax->f(ii,jj)=wcoax->f(ii,jj)*multiplier;
			wmb->f(ii,jj)=wmb->f(ii,jj)*multiplier;
			wmbl->f(ii,jj)=wmbl->f(ii,jj)*multiplier;
			if (jj<=number) wca[ii][jj]=wca[ii][jj]*multiplier;
			if (ii==1 && jj<=number) {
				//rescale w5
				w5[jj]=w5[jj]*pow(rescalefactor,(PFPRECISION) jj);

				if (jj==number) {
					//rescale w3
					for (index=1;index<=number;index++) w3[index]=w3[index]*pow(rescalefactor,(PFPRECISION) (number-index+1));

				}
			}
		}
	}

	//rescale curE and prevE
	for (ii=((h<=(number-2))?1:(2*number-h-1));ii<=((h<=(number-2))?(number-h):number);ii++)
	for (dp=1;dp<=d-1;dp++){
		if (ii<number) {
			curE[dp][ii] = curE[dp][ii]*pow(rescalefactor,(PFPRECISION)(dp+1));
			prevE[dp][ii+1] = prevE[dp][ii+1]*pow(rescalefactor,(PFPRECISION)(dp+1));
		}
	}

	//rescale datatable
	data->rescaledatatable(rescalefactor);
}


//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w3
//void rescaleatw3(int ii, structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
//				 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {



//}

//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w5
//void rescaleatw5(int jj,structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
//				 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {

	//rescale the previously filled arrays
//	rescale(1, jj, ct, data, v, w, wl, wcoax, wmb, wmbl, w5, w3, rescalefactor);

//	w5[jj]=w5[jj]*pow(rescalefactor,(double) jj);

//}

//rescale the entries in datatable
void pfdatatable::rescaledatatable(PFPRECISION rescalefactor) {

	scaling=scaling*rescalefactor;
	int i,j,k,l,m,n,o,p;


	for (i=0;i<31;i++) {
		inter[i] = inter[i]*pow(rescalefactor,i+2);
		bulge[i] = bulge[i]*pow(rescalefactor,i+2);
		hairpin[i] = hairpin[i]*pow(rescalefactor,i+2);

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					dangle[i][j][k][l] = dangle[i][j][k][l]*rescalefactor;
				}
				for (l=0;l<6;l++) {
					stack[i][j][k][l]=stack[i][j][k][l]*pow(rescalefactor,2);

					tstackcoax[i][j][k][l]=tstackcoax[i][j][k][l]*pow(rescalefactor,2);

					tstack[i][j][k][l]=tstack[i][j][k][l]*pow(rescalefactor,2);
					tstkm[i][j][k][l]=tstkm[i][j][k][l]*pow(rescalefactor,2);

					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							iloop11[i][j][k][l][m][n]=iloop11[i][j][k][l][m][n]*pow(rescalefactor,4);
							for (o=0;o<6;o++) {
								iloop21[i][j][k][l][m][n][o]=iloop21[i][j][k][l][m][n][o]*pow(rescalefactor,5);
								for (p=0;p<6;p++) {
									iloop22[i][j][k][l][m][n][o][p]=iloop22[i][j][k][l][m][n][o][p]*pow(rescalefactor,6);
								}
							}


						}
					}
				}
			}
		}
	}

	for (i=0;i<=numoftloops;i++) {

		tloop[i] = tloop[i]*pow(rescalefactor,6);

	}

	for (i=0;i<=numoftriloops;i++) {

		triloop[i] = triloop[i]*pow(rescalefactor,5);
	}

	for (i=0;i<=numofhexaloops;i++) {

		hexaloop[i] = hexaloop[i]*pow(rescalefactor,8);

	}



}

//thresh-structure builds a structure containing all base pairs above the probability thresh (expressed as a fraction from 0 to 1).
//Note that thresh must be 0.5 or larger for the resulting structure to be a valid secondary structure.
void thresh_structure(structure *ct, char *pfsfile, double thresh) {


	int i,j;
	short vers;
	PFPRECISION *w5, *w3, scaling;
	pfunctionclass *v, *w,*wmb,*wmbl,*wl,*wcoax;
	forceclass *fce;
	bool *mod,*lfce;
	pfdatatable *data;

	//allocate the ct file by reading the save file:
	ifstream sav(pfsfile,ios::binary);


	read(&sav,&(vers));//read the version of the save file
		//right now there is no infrastructure to indicate the wrong version is being read.
		//This should be changed in the future...

	read(&sav,&(ct->numofbases));

	sav.close();

	//allocate everything

	//array = new PFPRECISION *[ct.numofbases+1];
	//for (i=1;i<=ct.numofbases;i++) {
	//	array[i] = new PFPRECISION [i+1];
	//}
	//for (i=1;i<=ct.numofbases;i++) {
	//	for (j=1;j<=i;j++) {
	//		array[i][j]=0;
	//	}
	//}

	ct->allocate(ct->numofbases);

	w = new pfunctionclass(ct->numofbases);
	v = new pfunctionclass(ct->numofbases);
	wmb = new pfunctionclass(ct->numofbases);
	wmbl = new pfunctionclass(ct->numofbases);
	wcoax = new pfunctionclass(ct->numofbases);
	wl = new pfunctionclass(ct->numofbases);
	fce = new forceclass(ct->numofbases);

	w5 = new PFPRECISION [ct->numofbases+1];
	w3 = new PFPRECISION [ct->numofbases+2];

	lfce = new bool [2*ct->numofbases+1];
    mod = new bool [2*ct->numofbases+1];

	data = new pfdatatable();

	//load all the data from the pfsavefile:
	readpfsave(pfsfile, ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&scaling,mod,lfce,data);

	//reset the base pairing info:
	for (i=1;i<=ct->numofbases;i++) ct->basepr[1][i]=0;


	//fill array with the values for the plot:
	for (i=1;i<ct->numofbases;i++) {
		for (j=i+1;j<=ct->numofbases;j++) {



			if(calculateprobability(i,j,v,w5,ct,data,lfce,mod,scaling,fce)>thresh) {

				ct->basepr[1][i]=j;
				ct->basepr[1][j]=i;

			}



		}
	}

	//now build the structure

	ct->numofstructures = 1;

	delete w;
	delete v;
	delete wmb;
	delete fce;
	delete[] w5;
	delete[] w3;
	//for (i=1;i<=ct.numofbases;i++) {
	//	delete[] array[i];
	//}
	//delete[] array;
	delete[] lfce;
	delete[] mod;
	delete data;
	delete wmbl;
	delete wl;
	delete wcoax;


}
