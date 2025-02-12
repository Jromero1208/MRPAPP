//-*-C++-*-

#ifndef PAIRING_H
#define PAIRING_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "rpa.h"
#include "rpa_CuO.h"
#include "sepBasis.h"
#include "utilities.h"
// #include "bands.h"
#include "bandstructure.h"
#include "chi0.h"
#include "interpolation.h"
#include "ferminator.h"



namespace rpa {
	extern "C" void 
#ifdef glyph
	dgeev_
#else
	dgeev
#endif
	// JOBVL,JOBVR, N,    A,        LDA, WR,      WI,      VL,       LDVL,   VR,      , LDVR,  WORK,   LWORK, INFO
	(char *,char *,int *, double *,int *,double *,double *, double *, int *, double *, int *, double *, int *, int *);
	inline void GEEV(char jobvl,char jobvr,int n,psimag::Matrix<double> &A, int lda,
					 std::vector<double>  &wr,std::vector<double>  &wi, psimag::Matrix<double> &vl,
					 int ldvl, psimag::Matrix<double> &vr, int ldvr, std::vector<double> &work, int &lwork, int *info) 
					 {
#ifdef glyph						
						dgeev_
#else
						dgeev
#endif
						(&jobvl,&jobvr,&n,&(A(0,0)),&lda,
						      &(wr[0]),&(wi[0]),&(vl(0,0)),
						      &ldvl, &(vr(0,0)),&ldvr,&(work[0]),&lwork,info);
						}


	template<typename Field, typename BandsType, 
	         typename SuscType, typename GapType, template<typename> class MatrixTemplate, 
	         typename ConcurrencyType>
	class pairing {


	private:
		typedef Field 							FieldType;
		typedef std::complex<Field>				ComplexType;
		typedef MatrixTemplate<Field> 			MatrixType;
		typedef MatrixTemplate<size_t> 			IntMatrixType;
		typedef MatrixTemplate<ComplexType> 	ComplexMatrixType;
		typedef std::vector<Field>      		VectorType;
		typedef std::vector<ComplexType>      	ComplexVectorType;

		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		const size_t& dim;
		ConcurrencyType& conc;
		size_t interpolateChi_;
		size_t storeChi_;
		size_t readChi_;
		momentumDomain<FieldType,psimag::Matrix,ConcurrencyType>& qMesh;
		const size_t& nOrb,nk;
		size_t msize;
		momentumDomain<FieldType,psimag::Matrix,ConcurrencyType> kMesh;
		BandsType bands;
		int nkFI;
		size_t iqcalc;
		std::vector<int> kFtoBand;
		std::vector<FieldType> deltakF;
		std::vector<FieldType> vkF;
		// std::vector<int> negKF;
		SuscType chi0;
		chi0q<FieldType,SuscType,BandsType,GapType,psimag::Matrix,ConcurrencyType> susq;

		// IntMatrixType kF1kF2ToQMinus;
		// IntMatrixType kF1kF2ToQPlus;

		// susceptibility<FieldType,psimag::Matrix,ConcurrencyType> chi0;
		interaction<FieldType,psimag::Matrix,ConcurrencyType> rpa;
		SuscType chiRPAs;
		SuscType chiRPAc;
		ComplexMatrixType temps;
		ComplexMatrixType tempc;
	
		MatrixType gammaPP;
		FieldType paritySign;
		VectorType normalization;

		size_t nTotal;
		std::vector<SuscType> chiStore;
		std::vector<VectorType> qStore;
		ferminator<FieldType,BandsType,MatrixTemplate,ConcurrencyType> FSpoints;
		size_t nkF;
		MatrixType chikk;


	public:

		pairing(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters, 
				ConcurrencyType& concurrency, 
				const size_t interpolateChi,
				momentumDomain<FieldType,psimag::Matrix,ConcurrencyType>& qMeshIn):
			param(parameters),
			dim(param.dimension),
			conc(concurrency),
			interpolateChi_(interpolateChi),
			storeChi_(param.storeChi),
			readChi_(param.readChi),
			qMesh(qMeshIn),
			nOrb(param.nOrb),
			nk(parameters.nkInt),
			msize(param.nOrb*param.nOrb),
			kMesh(param,conc,param.nkInt,param.nkIntz,param.dimension),
			bands(param,conc,kMesh,false),
			chi0(param,conc),
			susq(param,qMesh,param.chifile,conc),
			rpa(param),
			chiRPAs(param,conc),
			chiRPAc(param,conc),
			temps(msize,msize),
			tempc(msize,msize),
			gammaPP(0,0),
			paritySign(param.pairingSpinParity?-1.0:1.0), // 1=triplet,0=singlet
			normalization(0,0),
			nTotal(0),
			chiStore(0,SuscType(param,conc)),
			qStore(0,VectorType(3)),
			FSpoints(param,conc),
			nkF(FSpoints.nTotal),
			chikk(nkF,nkF)
		{			
			// determineKF(file);

			// size_t nkF(FSpoints.nTotal);
			nTotal = (nkF * nkF + nkF) / 2;
			// nTotal = nkF * nkF;
			gammaPP.resize(nkF,nkF);
			normalization.resize(nkF);
			// if (param.Case=="Emery") {
			// 	std::cout << "Setting parameters to zero!!!!!!!!!\n";
			// 	interpolateChi_ = 0;
			// 	storeChi_ = 0 ;
			// 	readChi_ = 0;
			// }
			if (interpolateChi_==0)	{
				kMesh.set_momenta(false); // for chi claculation
				if (storeChi_==1 || readChi_==1) {
					chiStore.resize(nTotal,SuscType(param,conc));
					qStore.resize(nTotal,VectorType(3));
				}
				if (readChi_ == 1) {
					if (conc.rank()==0) readChiqTxt2(qStore,chiStore,param.chifile);
					std::cout << "ChiFile was read in \n";
					for (size_t iq=0;iq<qStore.size();iq++) {
						conc.broadcast(qStore[iq]);
						conc.broadcast(chiStore[iq]);
					}
					std::cout << "... and broadcasted \n";
					std::cout << "q,chi " << qStore[nTotal-1] << ", " << chiStore[nTotal-1].calcSus()<< "\n";
				}
			}
			calcNorm();
			calcGammaPP();
			if (interpolateChi_==0 && storeChi_==1 && readChi_==0) {
				for (size_t iq=0; iq<nTotal; iq++) conc.reduce(qStore[iq]);
				for (size_t iq=0; iq<nTotal; iq++) chiStore[iq].allReduce();
			}
			if (conc.rank()==0) {
				writeGammaPPAndNorm();
				calcEigenVectors();
				writeMatrixElementsOnFs();
				if (interpolateChi_==0 && storeChi_==1 && readChi_==0) {
					writeChiqTxt2(qStore,chiStore,"chiTest.txt");
				}
			}
		}


		void calcGammaPP() {
			// We have (nkF^2+nkF)/2 total elements to calculate for Gamma
			// size_t nTotal((nkF*nkF+nkF)/2);
			// size_t nkF(FSpoints.nTotal);
			std::vector<size_t> indToi(nTotal);
			std::vector<size_t> indToj(nTotal);
			size_t ind(0);
			for (size_t i = 0; i < nkF; ++i) for (size_t j = i; j < nkF; ++j){
			// for (size_t i = 0; i < nkF; ++i) for (size_t j = 0; j < nkF; ++j){
				indToi[ind] = i;
				indToj[ind] = j;
				ind++;
			}
			typedef PsimagLite::Range<ConcurrencyType> RangeType;
			RangeType range(0,nTotal,conc);
			VectorType Container(nTotal);
			VectorType Container2(nTotal);

			std::vector<FieldType> k1(3);
			std::vector<FieldType> k2(3);
			std::vector<FieldType> q(3);
			FieldType term1(0.0);
			FieldType chiTerm(0.0);
			// FieldType term2(0.0);

			for (;!range.end();range.next()) {
				ind = range.index();
				size_t ik1 = indToi[ind];
				size_t ik2 = indToj[ind];

				// if (ik1==0 && ik2==16) {

				k1[0] = FSpoints.kFx[ik1]; k1[1] = FSpoints.kFy[ik1]; k1[2] = FSpoints.kFz[ik1];
				k2[0] = FSpoints.kFx[ik2]; k2[1] = FSpoints.kFy[ik2]; k2[2] = FSpoints.kFz[ik2];
				size_t band1(FSpoints.kFtoBand[ik1]);
				size_t band2(FSpoints.kFtoBand[ik2]); 
				// q = kF1-kF2
				for (int l=0;l<3;l++) q[l] = k1[l]-k2[l];
				// for (int l=0;l<3;l++) q[l] = k2[l]-k1[l];
				if (param.Case == "Emery") calcGammaPPEmery(q,k1,k2,ik1,ik2,band1,band2,term1);
				else calcGammaPPTerms(ind,q,k1,k2,ik1,ik2,band1,band2,term1,chiTerm);
				Container[ind] = term1;
				Container2[ind] = chiTerm;
				if (conc.rank()==0) std::cout << "now calculating " << ind 
					                          << " of " << nTotal << 
					                             " terms with ik1=" << ik1 
					                          << " and ik2="<<ik2<<
					                             " term1=" << term1 << " , " 
					                          << "chiTerm=" << chiTerm <<
					                             // "; k1,k2,k2m=" << k1 << " ; " <<
					                             // k2 <<  " ; " << k2m <<
					                             "\n";
                 // }
			}

			conc.reduce(Container);
			conc.reduce(Container2);
			if (conc.rank()==0) {
				for (size_t ind = 0; ind < nTotal; ++ind){
					size_t i = indToi[ind];
					size_t j = indToj[ind];
					gammaPP(i,j) = Container[ind];
					gammaPP(j,i) = gammaPP(i,j);
					chikk(i,j) = Container2[ind];
					chikk(j,i) = chikk(i,j);
				}
				// MatrixType gammaPPFull(gammaPP);
				// for(size_t i=0;i<nkF;i++) for(size_t j=0;j<nkF;j++) {
				// 	gammaPP(i,j) = -0.5*(gammaPPFull(i,j)+paritySign*gammaPPFull(i,negKF[j]));
				// }
			}
		}

		void calcGammaPPEmery(std::vector<FieldType> q, 
							  VectorType& k1, VectorType& k2, size_t ik1,size_t ik2,
							  size_t band1, size_t band2, 
							  FieldType& result) {
			ComplexMatrixType chi0_gg(19,19);
			ComplexMatrixType chi0(3,3); //not needed (dummy)
			ComplexMatrixType chi0_g(19,3); //not needed (dummy)
			
			
			calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType> 
     	           calcChi0(param,kMesh,bands,q,conc,chi0,chi0_g,chi0_gg);

     	     // Now calculate RPA result of interaction
			interactionEmery<FieldType,psimag::Matrix,ConcurrencyType> rpaEmery(param);

				// First calculate the effective interaction with chi0
			ComplexMatrixType GammaS(19,19);
			ComplexMatrixType GammaC(19,19);
			ComplexMatrixType bareSpin(19,19);
			ComplexMatrixType bareCharge(19,19);

			rpaEmery.calcRPAResult(chi0_gg,1,GammaS,bareSpin,q); // renormalized spin interaction
			
			rpaEmery.calcRPAResult(chi0_gg,0,GammaC,bareCharge,q); // renormalized charge interaction


			// Determine GammaC and GammaS from chiRPA
			ComplexMatrixType chiRPAC(19,19);
			ComplexMatrixType chiRPAS(19,19);
			rpaEmery.calcChiRPAFromGammaRPA(GammaC,0,q,chiRPAC);
			rpaEmery.calcChiRPAFromGammaRPA(GammaS,1,q,chiRPAS);
			// Now re-calculate GammaC from chiRPA
			rpaEmery.setupVBare(q);
			rpaEmery.calcGammaFromChiRPA(rpaEmery.V_Charge,rpaEmery.V_Charge_coupl, chiRPAC, GammaC);
			rpaEmery.calcGammaFromChiRPA(rpaEmery.V_Spin  ,rpaEmery.V_Spin_coupl  , chiRPAS, GammaS);
			

			calcGammaPPOrbEmery(GammaS,GammaC,bareSpin,bareCharge,k1,k2,temps);

			calcGammaPPBand(temps,k1,k2,band1,band2,result);
		}

		void calcGammaPPOrbEmery(ComplexMatrixType& GammaS, ComplexMatrixType& GammaC, 
			                     ComplexMatrixType& bareSpin, ComplexMatrixType& bareCharge, 
								 VectorType& k1, VectorType& k2, 
								 ComplexMatrixType& result) {
			// the singlet interaction is 1/2*(3*GammaS - GammaC)
			// which then is multiplied with the basis g^i_{l1,l2}(k) * g^j_{l3,l4}(-k')
			VectorType mk2(3,0);
			VectorType mk1(3,0);
			// FieldType pairingFromCharge(1);
			// FieldType pairingFromSpin(0);
			for (size_t i=0;i<mk2.size();++i) mk2[i] = -k2[i];	
			for (size_t i=0;i<mk1.size();++i) mk1[i] = -k1[i];	
			sepBasis<FieldType,psimag::Matrix,ConcurrencyType> basisLeft(param,conc,k2);
			sepBasis<FieldType,psimag::Matrix,ConcurrencyType> basisRight(param,conc,mk1);
			// sepBasis<FieldType,psimag::Matrix,ConcurrencyType> basisLeft(param,conc,mk1);
			// sepBasis<FieldType,psimag::Matrix,ConcurrencyType> basisRight(param,conc,k2);

			for (size_t l1=0; l1<nOrb; l1++) for (size_t l2=0; l2<nOrb; l2++) {
			for (size_t l3=0; l3<nOrb; l3++) for (size_t l4=0; l4<nOrb; l4++) {
				size_t ind1 = l2+l1*nOrb;	
				size_t ind2 = l4+l3*nOrb;	
				result(ind1,ind2) = ComplexType(0.0,0.0);
				for (size_t i=0;i<19;++i) for (size_t j=0;j<19;++j) {
					result(ind1,ind2) += basisLeft(i,l1,l2) * basisRight(j,l3,l4) * 0.5 *
					// result(ind1,ind2) += basisLeft(i,l3,l4) * basisRight(j,l1,l2) *
					                    (
					                      (- 3 * (GammaS(i,j)-bareSpin(i,j))   - param.staticUFactor*bareSpin(i,j) )  * FieldType(param.pairingFromSpin) +
					                      (      (GammaC(i,j)-bareCharge(i,j)) + param.staticUFactor*bareCharge(i,j)) * FieldType(param.pairingFromCharge)
   					                    ); 
					                     // 1.0*(GammaC(i,j)); // only consider charge channel
				}
			}
			}
			return;
		}

		void calcGammaPPTerms(size_t ind,std::vector<FieldType> q, 
							  VectorType& k1, VectorType& k2, size_t ik1,size_t ik2,
							  size_t band1, size_t band2, 
							  FieldType& result,
							  FieldType& chiTerm) {
			// SuscType* chiq;
			SuscType chiq(param,conc);
			SuscType chi0Intq(param,conc);
			if (interpolateChi_==1) {
				getChi0forQ(q,chi0Intq);
				chiq = chi0Intq;
			} else if (interpolateChi_==0 && readChi_==0) {
				// std::cout << "Calculating chi for q=" << q << "\n";
				
				calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType> 
				          calcChi0(param,kMesh,bands,q,conc,chi0);
				// for (size_t i=0;i<msize;i++) for (size_t j=0;j<msize;j++) chi0(i,j) = calcChi0(i,j);
				// chi0.setLowerTriangle();
				chiq = chi0;
				if (storeChi_==1) {
					qStore[ind] = q;
					chiStore[ind] = chi0;
				}
			} else if (interpolateChi_==0 && readChi_==1) {
					chiq = chiStore[ind];
			}
			
			rpa.calcRPAResult(chiq,rpa.spinMatrix,chiRPAs,q);
			chiTerm = real(chiRPAs.calcSus());
			// if (ik1==0 && ik2==24) std::cout << "in calcGammaPPTerms: chiRPAs=" << chiRPAs.calcSus() << "\n";
			rpa.calcRPAResult(chiq,rpa.chargeMatrix,chiRPAc,q);
			// if (conc.rank()==0) std::cout << "in calcGammaPPTerms: chiRPAc=" << chiRPAc.calcSus() << "\n";
			matMul(chiRPAs,rpa.spinMatrix,temps);
			matMul(chiRPAc,rpa.chargeMatrix,tempc);
			matMul(rpa.spinMatrix,temps,chiRPAs);
			matMul(rpa.chargeMatrix,tempc,chiRPAc);
			// std::vector<std::complex<double> > chiRow(25,0);
			// if (ik1==0&&ik2==24) std::cout << "in calcGammaPPTerms chiRPAs=" << chiRPAs.calcSus() << "\n";
			calcGammaPPOrb(rpa.spinMatrix,chiRPAs,rpa.chargeMatrix,chiRPAc,temps);

			calcGammaPPBand(temps,k1,k2,band1,band2,result);
			// if (ik1==0&&ik2==40) std::cout << "in calcGammaPPTerms result=" << result << "\n";
		}

		void calcGammaPPOrb(const ComplexMatrixType& ms,
     						const ComplexMatrixType& ucsu, 
     						const ComplexMatrixType& mc,
     						const ComplexMatrixType& uccu,
     						      ComplexMatrixType& result) {
										  	
			for (size_t l1=0; l1<nOrb; l1++) {
				for (size_t l2=0; l2<nOrb; l2++) {
					for (size_t l3=0; l3<nOrb; l3++) {
						for (size_t l4=0; l4<nOrb; l4++) {
							size_t ind1 = l2+l1*nOrb;	
							size_t ind2 = l4+l3*nOrb;	
							size_t ind3 = ind2;
							size_t ind4 = ind1;
							if (param.pairingSpinParity==0) { // Singlet vertex
								result(ind1,ind2) = 0.5*(ms(ind3,ind4)-mc(ind3,ind4))*param.staticUFactor 
							                      - 0.5*uccu(ind3,ind4)*param.chargeFactor
							                      + 3./2.*ucsu(ind3,ind4);
		                    } else if (param.pairingSpinParity==1) { // Triplet vertex
								result(ind1,ind2) =  -0.5*uccu(ind3,ind4) - 0.5*ucsu(ind3,ind4);
							}
						}
					}
				}
			}
		}

		void calcGammaPPBand(const ComplexMatrixType& gammaOrb, 
							 VectorType& k1, 
							 VectorType& k2, 
							 size_t band1,
							 size_t band2,
							 FieldType& result) {

			VectorType ek1(nOrb,0);
			VectorType ek2(nOrb,0);
			ComplexMatrixType ak1(nOrb,nOrb);
			ComplexMatrixType ak2(nOrb,nOrb);
			ComplexMatrixType ak1m(nOrb,nOrb);
			ComplexMatrixType ak2m(nOrb,nOrb);

			VectorType k(3);
			VectorType km(3);
			for (size_t i=0; i<3; i++) km[i]=-k1[i];
			bands.getEkAndAk(k1,ek1,ak1);
			bands.getEkAndAk(km,ek1,ak1m,-1);
			for (size_t i=0; i<3; i++) km[i]=-k2[i];
			bands.getEkAndAk(k2,ek2,ak2);
			bands.getEkAndAk(km,ek2,ak2m,-1);


			result = FieldType(0.0);
			ComplexType c2(0.0,0.0);
			for (size_t l1=0; l1<nOrb; l1++) for (size_t l2=0; l2<nOrb; l2++) {
			for (size_t l3=0; l3<nOrb; l3++) for (size_t l4=0; l4<nOrb; l4++) {
				size_t ind1 = l2+l1*nOrb;	
				size_t ind2 = l4+l3*nOrb;	


				// std::cout << "l1,l2,l3,l4,band1,band2,aks: " << l1<<","<<l2<<","<<l3<<","<<l4<<","<<
				              // band1<<","<<band2<<","<<ak1(l2,band1)<<" , "<<ak1m(l3,band1)<<" , " << ak2(l1,band2) <<" , "<< ak2m(l4,band2) << "\n";

        		// ComplexType c1 = conj(ak1(l2,band1)) * conj(ak1m(l3,band1)) * 
		              			 // ak2(l1,band2) * ak2m(l4,band2);
        		
        		// ComplexType c1 = conj(ak1(l2,band1)) * ak1(l3,band1) * 
		              			 // ak2(l1,band2) * conj(ak2(l4,band2));
				
				ComplexType c1 = ak2(l2,band2) * ak2m(l3,band2)*     // works for Emery model and general case!!!
                                 conj(ak1(l1,band1)) * conj(ak1m(l4,band1)); // as in Kreisel et al. PRB 88, 094522

                c2 += c1*gammaOrb(ind1,ind2);
                
                // c2 += conj(c1)*gammaOrb(ind1,ind2);

                // if (ind1==8 && ind2==4) std::cout << "c1,k1,k2,a...a" << c1<<ak2(l2,band2) <<  ak2m(l3,band2)
                                 // << conj(ak1(l1,band1)) << conj(ak1m(l4,band1)) << "\n";
		        // result += real(c1*gammaOrb(ind1,ind2));
			}
			}
			result = real(c2);
			// if (conc.rank()==0) std::cout << "c2 = " << c2 << "\n";
						 	
		 }

		void getChi0forQ(VectorType& q, SuscType& chi) {
			// mapQto1Quadrant(q);
			// ComplexVectorType susVec(qMesh.nktot);
			interpolation<FieldType,MatrixTemplate,SuscType,ConcurrencyType> interpolate(qMesh,susq);
			if (param.dimension==2) {
				interpolate.BiLinearGeneral(q,chi);
				// interpolate.BiLinear(q,chi);
			}
			else if (param.dimension==3) {
				interpolate.TriLinearGeneral(q,chi);
			}
		}



		void calcNorm() {
			size_t nkF(FSpoints.nTotal);
			for (size_t ik=0; ik<nkF; ++ik) {
				VectorType k(3,0);
				// k[0] = kFx[ik]; k[1]=kFy[ik]; k[2]=kFz[ik];
				// size_t band(kFtoBand[ik]);
				// FieldType vkF2 = bands.fermiVelocity2D(k,band);
				// std::cout << k[0] << " , " << k[1] << " , " << band << " , " << vkF2 << "\n";
				normalization[ik]=FSpoints.deltakF[ik]/(pow(2.*param.pi_f,dim)*FSpoints.vkF[ik]);
			}
		}

		void writeGammaPPAndNorm() {
			// std::cout << "kFx1,kFy1" << kFx[0] << "," << kFy[0];
			// First in a format suitable for reading into R using the script in gammapp_diag.R
			size_t nkF(FSpoints.nTotal);
			std::ofstream os("Gammakkp.txt");
			os << "nkF:\n" << nkF << "\n";
			os << "kFx:\n" << FSpoints.kFx << "\n";
			os << "kFy:\n" << FSpoints.kFy << "\n";
			os << "kFz:\n" << FSpoints.kFz << "\n";
			os << "U,U',J,J':\n" << param.U << " , " << param.Up << " , " << param.J << " , " << param.Jp << "\n\n";
			os << "Gamma(k,k'): \n";
			os << gammaPP << "\n\n";
			os << "Normalization: \n";
			os << normalization << "\n\n";
			os.close();

			// Now write in jsn format for post-processing with python/matplotlib
			std::ofstream os2("Gammakkp.jsn");
			int width(13);
			os2.precision(width);
			os2 << std::fixed; // scientific;
			os2 << "{ ";
			os2 << "\"U\": "  << param.U  << ",\n";
			os2 << "\"Up\": " << param.Up << ",\n";
			os2 << "\"J\": "  << param.J  << ",\n";
			os2 << "\"Jp\": " << param.Jp << ",\n";
			os2 << " \"numberKfPoints\": " << nkF << ",\n";
			os2 << " \"kfPoints\": [\n";
			for (size_t ik=0;ik<nkF;ik++) {
				os2 << "[" << FSpoints.kFx[ik] << " ," << FSpoints.kFy[ik] << " ," << FSpoints.kFz[ik] << "] ";
				if (ik<nkF-1) os2 << ",";
			}
			os2 << "],\n";
			os2 << " \"GammaPP\": [\n";
			for(size_t ik=0;ik<nkF;ik++) {
				os2 << "[";
				for(size_t ikp=0;ikp<nkF;ikp++) {
					os2 << gammaPP(ik,ikp);
					if (ikp<nkF-1) os2 << ", ";
				}
				os2 << "]";
				if (ik<nkF-1) os2 << ",\n";
			}
			os2 << "],\n";
			os2 << " \"Measure\": [\n";
			os2 << "[";
			for(size_t ik=0;ik<nkF;ik++) {
				os2 << normalization[ik];
				if (ik<nkF-1) os2 << ", ";
			}
			os2 << "]],\n";
			os2 << " \"chikk\": [\n";
			for(size_t ik=0;ik<nkF;ik++) {
				os2 << "[";
				for(size_t ikp=0;ikp<nkF;ikp++) {
					os2 << chikk(ik,ikp);
					if (ikp<nkF-1) os2 << ", ";
				}
				os2 << "]";
				if (ik<nkF-1) os2 << ",\n";
			}
			os2 << "]\n";
			os2 << "}\n";
			os2.close();

		}

	void eigen(MatrixType& matrix, VectorType& wr, MatrixType& vr) const {
		int n = matrix.n_row();
		int lwork = 4*n;
		std::vector<double> work(lwork);
		int info;
		std::vector<double> wi(n);
		psimag::Matrix<double> vl(n,n);


		GEEV('N','V',n,matrix,n,wr,wi,vl,n,vr,n,work,lwork,&info);
		if (info!=0) {
			throw std::runtime_error("dgeev: failed\n");
		}
	}

	void calcEigenVectors() {
		std::cout << "Now calculating the eigenvalues and -vectors of BSE\n";
		MatrixType matrix(gammaPP);
		MatrixType eigenvects(gammaPP);
		for (size_t ik1=0;ik1<nkF;ik1++) {
			for (size_t ik2=0;ik2<nkF;ik2++) {
				matrix(ik1,ik2) = -gammaPP(ik1,ik2) * normalization[ik2];
			}
		}
		VectorType eigenvals(nkF);
		eigen(matrix,eigenvals,eigenvects);
		// Now print out eigenvalues and -vectors
		std::ofstream os2("Gap.jsn");
		int width(13);
		os2.precision(width);
		os2 << std::fixed; // scientific;
		os2 << "{ ";
		os2 << "\"U\": "  << param.U  << ",\n";
		os2 << "\"Up\": " << param.Up << ",\n";
		os2 << "\"J\": "  << param.J  << ",\n";
		os2 << "\"Jp\": " << param.Jp << ",\n";
		os2 << " \"numberKfPoints\": " << nkF << ",\n";

		os2 << " \"kfPoints\": [\n";
		for (size_t ik=0;ik<nkF;ik++) {
			os2 << "[" << FSpoints.kFx[ik] << " ," << FSpoints.kFy[ik] << " ," << FSpoints.kFz[ik] << "] ";
			if (ik<nkF-1) os2 << ",";
		}
		os2 << "],\n";

		os2 << " \"Eigenvectors\": [\n";
		for(size_t ikp=0;ikp<nkF;ikp++) {
			os2 << "[";
			for(size_t ik=0;ik<nkF;ik++) {
				os2 << eigenvects(ik,ikp);
				if (ik<nkF-1) os2 << ", ";
			}
			os2 << "]\n";
			if (ikp<nkF-1) os2 << ",\n";
		}

		os2 << "],\n";
		os2 << " \"Eigenvalues\": [\n";
		os2 << "[";
		for(size_t ik=0;ik<nkF;ik++) {
			os2 << eigenvals[ik];
			if (ik < nkF-1) os2 << ", ";
		}
		os2 << "]]\n";
		os2 << "}\n";
		os2.close();

	}

	void writeMatrixElementsOnFs(){
		std::ofstream os("akOnFS.txt");
		os << "nkF: \n";
		os << nkF << "\n";

		for (size_t ik=0;ik<nkF;ik++) {
			size_t band(FSpoints.kFtoBand[ik]);
			VectorType k(3,0);
			k[0] = FSpoints.kFx[ik];
			k[1] = FSpoints.kFy[ik];
			k[2] = FSpoints.kFz[ik];

			ComplexMatrixType ak(nOrb,nOrb);
			VectorType ek(nOrb);
			bands.getEkAndAk(k,ek,ak);

			os << k[0] << " " << k[1] << " " << k[2] << " ";
			for (size_t l1=0;l1<nOrb;l1++) for (size_t l2=0;l2<nOrb;l2++) {
				ComplexType c1(ak(l1,band)*conj(ak(l2,band)));
				if (imag(c1)>=0) os << real(c1) << "+" << imag(c1) << "j" << " ";
				else  os << real(c1) << imag(c1) << "j" << " ";
			}
			os << "\n";
		}

	}

	};


}

#endif
