//-*-C++-*-

#ifndef SUSCEPTIBILITY_H
#define SUSCEPTIBILITY_H

#include <string>
#include <vector>
#include <fstream>
// #include "math.h"
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "Fermi.h"
#include "rpa.h"
#include "bandstructure.h"
#include "chi0.h"
#include "gaps2D.h"
#include "gaps3D.h"
#include "ferminator.h"
#include "rpa_CuO.h"
#include "greensFunction.h"


namespace rpa {

	template<typename Field, typename SuscType, typename BandsType,
	         template<typename> class MatrixTemplate,
	         typename ConcurrencyType>
	class susceptibility {
		private:
			typedef Field 				                     FieldType;
			typedef std::complex<Field>		             ComplexType;
			typedef MatrixTemplate<Field> 			       MatrixType;
			typedef MatrixTemplate<ComplexType>    	   ComplexMatrixType;
			typedef std::vector<Field>      	    	   VectorType;
			typedef std::vector<std::complex<Field> >  ComplexVectorType;
			typedef PsimagLite::Range<ConcurrencyType> RangeType;
			typedef std::vector<SuscType>              VectorSuscType;
#ifdef USE_SCGAP2D
			typedef rpa::gap2D<FieldType,psimag::Matrix,ConcurrencyType> GapType;
#elif  USE_SCGAP3D
			// typedef rpa::gap3D<FieldType,psimag::Matrix,BandsType,ConcurrencyType> GapType;
			typedef rpa::gap3D<FieldType,psimag::Matrix,ConcurrencyType> GapType;
#else
			typedef rpa::gap2D<FieldType,psimag::Matrix,ConcurrencyType> GapType;
#endif

			const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
			ConcurrencyType& conc;
			// momentumDomain<Field,psimag::Matrix> qMesh;
			size_t numberOfQ;
			size_t msize;
			size_t nq1;
			size_t nq2;
			size_t nq3;
			size_t nw;
			VectorType omega;
			std::vector<std::vector<FieldType> > QVec;
			std::vector<size_t> indexToiq;
			std::vector<size_t> indexToiw;
			FieldType wmin_,wmax_;
			std::vector<size_t> qtok;

			greensFunction<FieldType,psimag::Matrix,ConcurrencyType> g;


		public:
			// typedef std::vector<SuscType> BaseType;

// Class Constructor
		susceptibility(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
					   ConcurrencyType& concurrency,
					   const FieldType& qxmin, const FieldType& qxmax, const size_t nq1In,
					   const FieldType& qymin, const FieldType& qymax, const size_t nq2In,
					   const FieldType& qzmin, const FieldType& qzmax, const size_t nq3In,
					   const FieldType& wmin,  const FieldType& wmax,  const size_t nwIn
					   ):
						param(parameters),
						conc(concurrency),
						// qMesh(param,nqxIn,nqzIn,3),
						numberOfQ(nq1In*nq2In*nq3In*nwIn),
						msize(param.nOrb*param.nOrb),
						nq1(nq1In),
						nq2(nq2In),
						nq3(nq3In),
						nw(nwIn),
						omega(nw,0),
						QVec(numberOfQ,VectorType(5,0)),
						indexToiq(numberOfQ,0),
						indexToiw(numberOfQ,0),
						wmin_(wmin),
						wmax_(wmax),
						qtok(numberOfQ,0),
						g(param,conc)
		{

			if (conc.rank()==0) std::cerr << "numberOfQ: " << numberOfQ << "\n";



				// Setup q-mesh based on nq's and q-limits
			if(param.Green) {
				g.init();
				if (param.Gfile == "") g.printGreens();
				if (param.Qmesh){
					setupQandOmegaMesh(nq1,nq2,nq3,numberOfQ,nw,wmin,wmax,QVec,qtok,g);
				} else {
					setupQandOmegaMesh(nq1,nq2,nq3,numberOfQ,nw,wmin,wmax,QVec,qtok,g,param.Qmesh);
				}
			} else {
				setupQandOmegaMesh(nq1,nq2,nq3,numberOfQ,nw,
				               qxmin,qxmax,qymin,qymax,qzmin,qzmax,
				               wmin,wmax,QVec);
			}

			if (param.scState==1 && param.printGap==1 && conc.rank()==0)
			{
				if (conc.rank()==0) std::cout << "Now writing gap.txt \n";
				printGap3();
				if (conc.rank()==0) std::cout << "Done writing gap.txt \n";
			}



			if (param.Case=="Emery") {
				std::vector<ComplexMatrixType> chi0(numberOfQ,ComplexMatrixType(3,3));
				std::vector<ComplexMatrixType> chi0_g(numberOfQ,ComplexMatrixType(19,3));
				std::vector<ComplexMatrixType> chi0_gg(numberOfQ,ComplexMatrixType(19,19));

				std::string filename("chi0Emery.txt");
				if (param.readChiForSus) readChi0Emery(QVec,chi0,chi0_g,chi0_gg,filename);
				else calcEmeryChi0(chi0,chi0_g,chi0_gg);
				calcEmeryRPA(chi0,chi0_g,chi0_gg);

			} else {

				std::cerr << "Now Calculating Chi0" << std::endl;
				VectorSuscType ChiMat(nq1In*nq2In*nq3In*nwIn,SuscType(parameters,concurrency));
				std::vector<std::vector<VectorSuscType> > chi0Matrix(param.nSite, std::vector<VectorSuscType> (param.nSite, ChiMat));
				/*Average Chi0 code

				for (size_t i = 0; i < param.NumTB; i++) {
					chi0Matrix.push_back(ChiMat);

					if (param.Green){
						calcElements(chi0Matrix[i],param.Green);
					} else {
						calcElements(chi0Matrix[i],param.tbfile[i]);
					}

					// Average Chiq
					for (size_t j = 0; j < nq1In*nq2In*nq3In*nwIn; j++) {
						for (size_t k = 0; k < (param.nOrb*param.nOrb); k++) {
							for (size_t l = 0; l < (param.nOrb*param.nOrb); l++) {
								ChiMat[j](k,l) += chi0Matrix[i][j](k,l)/(double)param.NumTB;
							}
						}
					}

				}*/
				// make the split for multisite calculation
				if (param.nSite == 1) {
					if (param.Green){
						calcElements(chi0Matrix[0][0],param.Green,param.tbfile);
					} else {
						calcElements(chi0Matrix[0][0],param.tbfile);
					}
					if (conc.rank()==0) std::cout << "Now printing out chiq \n"; writeChiqTxt(chi0Matrix[0][0]);
				} else {
					calcElements(chi0Matrix,param.tbfile,param.nSite);


					//print out Chiq and chiRPA for optimized system
					if (conc.rank()==0) std::cout << "Now printing out chiq \n"; writeChiqTxt(chi0Matrix,param.nSite);

				}


			}

		}

		//standard calcElements

		void calcElements(VectorSuscType& chi0Matrix, std::string file) {

				// Setup k-mesh for chi0 calculation
			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			kmesh.set_momenta(false);
			BandsType bands(param,conc,kmesh,false,file);



			RangeType range(0,QVec.size(),conc);

			// SuscType chi0QW(param,conc);
			for (;!range.end();range.next()) {

				size_t iQ = range.index();
				std::vector<FieldType> q(3);
				q[0]=QVec[iQ][0]; q[1]=QVec[iQ][1]; q[2]=QVec[iQ][2];

				if (param.scState==1)
				{
					GapType Delta(param,conc);
					calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
				               calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],Delta,QVec[iQ][3],0);
	      /*} else if (param.Green) {
					  calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
							//calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],iQ);
								calcChi0(param,bands,q,numberOfQ,conc,chi0Matrix[iQ],iQ,qtok[iQ],g);

				*/} else {
	          if (wmin_==0.0 && wmax_ == 0.0)
						{
	           				calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
				               calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],false,param.calcOnlyDiagonal);
	           } else {
		           			calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
					             calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],QVec[iQ][3]);
					   }
	      }

				if (conc.rank()==0) {
					std::cout.precision(7);
					std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3]
					          << "  of " << numberOfQ
	                          << " total. ChiPhys=" << chi0Matrix[iQ].calcSus()
	                          << "\n";
                }
			}

			for (size_t iq=0;iq<numberOfQ;iq++) chi0Matrix[iq].allReduce();
		}

		//standard calcElements for optimized system

		void calcElements(std::vector<std::vector<VectorSuscType> >& chi0Matrix, std::string file, int nSites) {

		//	std::cerr << "Calculating optimized Chi0" << std::endl;
				// Setup k-mesh for chi0 calculation
			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			kmesh.set_momenta(false);
			BandsType bands(param,conc,kmesh,false,file);

			std::vector< std::vector<int> > v;
			int vindex=0;
			for (size_t i = 0; i < param.nSite; i++) {
				std::vector<int> row;
				row.push_back(i);
				row.push_back(i);
				for (size_t j = i; j < param.nSite; j++) {
					row[1]=j;
					v.push_back(row);
					vindex++;
				}
			}

			RangeType range(0,QVec.size()*vindex,conc);
			// SuscType chi0QW(param,conc);
			for (;!range.end();range.next()) {

				size_t rindex = range.index();
				size_t iQ = rindex/vindex;
				std::vector<FieldType> q(3);
				q[0]=QVec[iQ][0]; q[1]=QVec[iQ][1]; q[2]=QVec[iQ][2];


					int i = v[rindex-iQ*vindex][0];
					int j = v[rindex-iQ*vindex][1];

					if (conc.rank() == 0) std::cerr << "now calculating x,y: " << i << " , " << j << std::endl;
					if (wmin_==0.0 && wmax_ == 0.0)
					{
						calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
						 		calcChi0(param,kmesh,bands,q,conc,chi0Matrix[i][j][iQ],i,j,false,param.calcOnlyDiagonal);
					 } else {
						 calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
								calcChi0(param,kmesh,bands,q,conc,chi0Matrix[i][j][iQ],i,j,QVec[iQ][3]);
					 }


						 if (conc.rank()==0) {
		 					std::cout.precision(7);
		 					std::cout << "iQ = " << iQ << " q= " << q
		 										<< "  of " << numberOfQ << " total. ChiPhys=" << chi0Matrix[i][j][iQ].calcSus(param.nOrbSite) << "\n";
		 								}


			}

			for (size_t iq=0;iq<numberOfQ;iq++){ for (size_t vin = 0; vin < vindex; vin++) {
				int i = v[vin][0];
				int j = v[vin][1];
				chi0Matrix[i][j][iq].allReduce();
			}
		}

		}

		//calcElements for green
		void calcElements(VectorSuscType& chi0Matrix, bool paramg, std::string file) {

			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			kmesh.set_momenta(false);
			BandsType bands(param,conc,kmesh,false,file);
			RangeType range(0,QVec.size(),conc);

			// SuscType chi0QW(param,conc);
			for (;!range.end();range.next()) {
				size_t iQ = range.index();
				std::vector<FieldType> q(3);
				q[0]=QVec[iQ][0]; q[1]=QVec[iQ][1]; q[2]=QVec[iQ][2];

				calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
				calcChi0(param,bands,q,numberOfQ,conc,chi0Matrix[iQ],iQ,qtok[iQ],g);

				calcChi0(param,kmesh,bands,q,conc,chi0Matrix[iQ],false,param.calcOnlyDiagonal);

				if (conc.rank()==0) {
					std::cout.precision(7);
					std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3]
					          << "  of " << numberOfQ << " total. ChiPhys=" << chi0Matrix[iQ].calcSus()
	                  << "\n";
        }
			}

			for (size_t iq=0;iq<numberOfQ;iq++) chi0Matrix[iq].allReduce();
		}




		void calcEmeryChi0(std::vector<ComplexMatrixType>& chi0,
						   std::vector<ComplexMatrixType>& chi0_g,
						   std::vector<ComplexMatrixType>& chi0_gg) {
			// Setup k-mesh for chi0 calculation
			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			kmesh.set_momenta(false);
			BandsType bands(param,conc,kmesh,false); // false = no Caching
			RangeType range(0,numberOfQ,conc);

			for (;!range.end();range.next()) {

				size_t iQ = range.index();
				std::vector<FieldType> q(3);
				q[0]=QVec[iQ][0]; q[1]=QVec[iQ][1]; q[2]=QVec[iQ][2];

   				calcChi0Matrix<FieldType,SuscType,BandsType,GapType,MatrixTemplate,ConcurrencyType>
				           calcChi0(param,kmesh,bands,q,conc,chi0[iQ],chi0_g[iQ],chi0_gg[iQ]);

				if (conc.rank()==0) std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3]
					          << "  of " << numberOfQ
	                          << "\n";
	           }


			for (size_t iq=0;iq<numberOfQ;iq++) {
				conc.allReduce(chi0[iq]);
				conc.allReduce(chi0_g[iq]);
				conc.allReduce(chi0_gg[iq]);
			}

			if (conc.rank()==0) {
				std::string	filename("chi0Emery.txt");
				writeChi0Emery(QVec,chi0,chi0_g,chi0_gg,filename);
			}
		}

		void calcEmeryRPA(std::vector<ComplexMatrixType>& chi0,
			              std::vector<ComplexMatrixType>& chi0_g,
			              std::vector<ComplexMatrixType>& chi0_gg) {
			// Now calculate the RPA Chi
			interactionEmery<FieldType,psimag::Matrix,ConcurrencyType> rpaEmery(param);
			std::vector<ComplexMatrixType> chiS(numberOfQ,ComplexMatrixType(3,3));
			std::vector<ComplexMatrixType> chiC(numberOfQ,ComplexMatrixType(3,3));

			for (size_t iQ=0;iQ<numberOfQ;iQ++) {
				// First calculate the effective interaction with chi0
				std::vector<FieldType> q(3);
				q[0]=QVec[iQ][0]; q[1]=QVec[iQ][1]; q[2]=QVec[iQ][2];
				ComplexMatrixType GammaS(19,19);
				ComplexMatrixType GammaC(19,19);
				ComplexMatrixType bareSpin(19,19);
				ComplexMatrixType bareCharge(19,19);
				ComplexMatrixType couplingSpin(19,19);
				ComplexMatrixType couplingCharge(19,19);
				rpaEmery.calcRPAResult(chi0_gg[iQ],1,GammaS,bareSpin,q); // renormalized spin interaction
				rpaEmery.calcRPAResult(chi0_gg[iQ],0,GammaC,bareCharge,q); // renormalized charge interaction
				// Now calculate the spin and charge susceptibilities
				ComplexMatrixType aGammaSa(3,3),aGammaCa(3,3);
				for (size_t l1=0;l1<3;l1++) for (size_t l2=0;l2<3;l2++) {
					for (size_t i=0;i<19;i++) for (size_t j=0;j<19;j++) {
						aGammaSa(l1,l2) += chi0_g[iQ](i,l1)*GammaS(i,j)*chi0_g[iQ](j,l2);
						aGammaCa(l1,l2) += chi0_g[iQ](i,l1)*GammaC(i,j)*chi0_g[iQ](j,l2);
					}
					chiS[iQ](l1,l2) = chi0[iQ](l1,l2) - aGammaSa(l1,l2);
					chiC[iQ](l1,l2) = chi0[iQ](l1,l2) - aGammaCa(l1,l2);
				}
				if (conc.rank()==0) {
					std::cout << "iQ = " << iQ << " q= " << q << " w = " << QVec[iQ][3];
					ComplexType susS = calcChiSRPA(chiS[iQ]);
					ComplexType susN = calcChiNematicRPA(chiC[iQ]);
					std::cout << "spin sus.: " << susS << " , nematic sus.: " << susN << "\n";
				}
			}
			writeChiSChiC(chiS,chiC);

		}

		ComplexType calcChiSRPA(ComplexMatrixType& chiS) {
			ComplexType sus(0.0);
			for (size_t i=0;i<chiS.n_row();i++) for (size_t j=0;j<chiS.n_row();j++) {
				sus += 0.5*chiS(i,j);
			}
			return sus;
		}

		ComplexType calcChiNematicRPA(ComplexMatrixType& chiC) {
			ComplexType sus(0.0);
			sus  = chiC(1,1)-chiC(1,2)+chiC(2,2)-chiC(2,1);
			return sus;
		}

		void setupQandOmegaMesh(size_t nq1, size_t nq2, size_t nq3,
								size_t numberOfQ,
								size_t nw,
								const FieldType& qxmin,
								const FieldType& qxmax,
								const FieldType& qymin,
								const FieldType& qymax,
								const FieldType& qzmin,
								const FieldType& qzmax,
								const FieldType& wmin,
								const FieldType& wmax,
								std::vector<std::vector<FieldType> >& QVec
								) {

			MatrixType momenta(nq1*nq2*nq3,3);
			if (nq1 == 1 && nq2 == 1 && nq3 == 1) { //Q is fixed
				momenta(0,0) = qxmin;
				momenta(0,1) = qymin;
				momenta(0,2) = qzmin;
			} else if (nq1 > 1 && nq2 == 1 && nq3 == 1) { // 1D Q-Scan
				for (size_t iq1=0;iq1<nq1;iq1++){
					momenta(iq1,0) = qxmin + float(iq1)/float(nq1-1) * (qxmax - qxmin);
					momenta(iq1,1) = qymin + float(iq1)/float(nq1-1) * (qymax - qymin);
					momenta(iq1,2) = qzmin + float(iq1)/float(nq1-1) * (qzmax - qzmin);
				}
			} else if (nq1 > 1 && nq2 > 1 && nq3 == 1) { // 2D Q-scan (in-plane)
				for (size_t iq2=0;iq2<nq2;iq2++){
					for (size_t iq1=0;iq1<nq1;iq1++){
						size_t index(iq1 + nq1*iq2);
						momenta(index,0) = qxmin + float(iq1)/float(nq1-1) * (qxmax - qxmin);
						momenta(index,1) = qymin + float(iq2)/float(nq2-1) * (qymax - qymin);
						momenta(index,2) = qzmin;
					}
				}
			} else if (nq1 > 1 && nq2 > 1 && nq3 > 1) { // 3D Q-scan
				for (size_t iq3=0;iq3<nq3;iq3++){
					for (size_t iq2=0;iq2<nq2;iq2++){
						for (size_t iq1=0;iq1<nq1;iq1++){
							size_t index(iq1 + nq1*iq2 + nq1*nq2*iq3);
							momenta(index,0) = qxmin + float(iq1)/float(nq1-1) * (qxmax - qxmin);
							momenta(index,1) = qymin + float(iq2)/float(nq2-1) * (qymax - qymin);
							momenta(index,2) = qzmin + float(iq3)/float(nq3-1) * (qzmax - qzmin);
						}
					}
				}
			}
			// Setup linear omega-mesh
			for (size_t i=0; i<nw; i++) omega[i] = wmin + (wmax - wmin) * float(i)/fmax(float(nw-1),float(1));
			// Now combine vectors
			for (size_t iq=0; iq<nq1*nq2*nq3; iq++) for (size_t iw=0; iw<nw; iw++) {
				size_t index = iw + iq * nw;
				indexToiq[index] = iq;
				indexToiw[index] = iw;
				QVec[index][0] = momenta(iq,0);
				QVec[index][1] = momenta(iq,1);
				QVec[index][2] = momenta(iq,2);
				QVec[index][3] = omega[iw];
			}


		}

		// for GreensFunction
		void setupQandOmegaMesh(size_t nq1, size_t nq2, size_t nq3,
								size_t numberOfQ,
								size_t nw,
								const FieldType& wmin,
								const FieldType& wmax,
								std::vector<std::vector<FieldType> >& QVec,
								std::vector<size_t>& qtok,
								greensFunction<FieldType,psimag::Matrix,ConcurrencyType> g)
								{

			MatrixType momenta(nq1*nq2*nq3,3);
			std::vector<size_t> qtokindex(numberOfQ,0);
			if (nq1 == 1 && nq2 == 1 && nq3 == 1) { //Q is fixed
				momenta(0,0) = g.momentumDomain1.momenta(0,0);
				momenta(0,1) = g.momentumDomain1.momenta(0,1);
				momenta(0,2) = g.momentumDomain1.momenta(0,2);
			} else {
				size_t nkz = g.momentumDomain1.nkz;
				size_t nk = g.momentumDomain1.nk;


				size_t index=0;
				for (size_t iq1=0;iq1<nq1;iq1++){
					size_t qtok1 = ((float)(iq1+nq1) / (float)(nq1+nq1)) * nk;


					for (size_t iq2=0;iq2<nq2;iq2++){
						size_t qtok2 = ((float)(iq2+nq2) / (float)(nq2+nq2)) * nk;

						for (size_t iq3=0;iq3<nq3;iq3++){
							size_t qtok3 = floor((float)(iq3+nq3) / (float)(nq3+nq3)) * nkz;
							//qtok3 = 0;
							index=(iq3 + iq2 * nq3 + iq1 * nq2 * nq3);

							qtokindex[index] = (qtok3 + qtok2*nkz + qtok1*nkz*nk);
							//std::cerr << "index, qtokindex[index] = " << index << ", " << qtokindex[index] << "," <<  (float)iq1/(float)nq1 << std::endl;
							momenta(index,0) = g.momentumDomain1.momenta(qtokindex[index],0);
							momenta(index,1) = g.momentumDomain1.momenta(qtokindex[index],1);
							momenta(index,2) = g.momentumDomain1.momenta(qtokindex[index],2);

							// shift index so vector addition is correct
							size_t qk1, qk2, qk3;
							if (qtok3 < (nkz/2)) {qk3 = qtok3 + (nkz/2);} else {qk3 = qtok3 - (nkz/2);}
							if (qtok2 < (nk/2)) {qk2 = qtok2 + (nk/2);} else {qk2 = qtok2 - (nk/2);}
							if (qtok1 < (nk/2)) {qk1 = qtok1 + (nk/2);} else {qk1 = qtok1 - (nk/2);}
							//qtok3 = 0;
							qtokindex[index] = (qk3 + qk2*nkz + qk1*nkz*nk);
							//index++;
						}
					}
				}
			}
			// Setup linear omega-mesh
			for (size_t i=0; i<nw; i++) omega[i] = wmin + (wmax - wmin) * float(i)/fmax(float(nw-1),float(1));
			// Now combine vectors
			for (size_t iq=0; iq<nq1*nq2*nq3; iq++) for (size_t iw=0; iw<nw; iw++) {
				size_t index = iw + iq * nw;
				indexToiq[index] = iq;
				indexToiw[index] = iw;
				QVec[index][0] = momenta(iq,0);
				QVec[index][1] = momenta(iq,1);
				QVec[index][2] = momenta(iq,2);
				QVec[index][3] = omega[iw];
				qtok[iq] = qtokindex[iq];
			}


		}

		// for GreensFunction if qmesh is off
		// sets a path of Q from 0, 0, 0, to Pi, Pi, Pi, or Pi, Pi, 0
		void setupQandOmegaMesh(size_t nq1, size_t nq2, size_t nq3,
								size_t numberOfQ,
								size_t nw,
								const FieldType& wmin,
								const FieldType& wmax,
								std::vector<std::vector<FieldType> >& QVec,
								std::vector<size_t>& qtok,
								greensFunction<FieldType,psimag::Matrix,ConcurrencyType> g,
								bool Qmesh)
								{

			MatrixType momenta(nq1*nq2*nq3,3);
			std::vector<size_t> qtokindex(numberOfQ,0);

				size_t nkz = g.momentumDomain1.nkz;
				size_t nk = g.momentumDomain1.nk;

				for (size_t iq2=0;iq2<= ((numberOfQ/2));iq2++){
					size_t iq3 = iq2 + (numberOfQ/2);

					size_t qtok3 = 0;
					if (nq3 > 1) qtok3 = ((float)iq3 / (float)nq3) * nkz;

						size_t qtok2 = ((float)iq3 / (float)numberOfQ) * nk;

							size_t qtok1 = ((float)iq3 / (float)numberOfQ) * nk;



							size_t index(iq2);
							qtokindex[index] = (qtok3 + qtok2*nkz + qtok1*nkz*nk);
							//std::cerr << "index, qtokindex[index] = " << index << ", " << qtokindex[index] << "," <<  (float)iq1/(float)nq1 << std::endl;
							momenta(index,0) = g.momentumDomain1.momenta(qtokindex[index],0);
							momenta(index,1) = g.momentumDomain1.momenta(qtokindex[index],1);
							momenta(index,2) = g.momentumDomain1.momenta(qtokindex[index],2);


							// shift index so vector addition is correct
							if (qtok3 < (nkz/2)) {qtok3 += (nkz/2);} else {qtok3 -= (nkz/2);}
							if (qtok2 < (nk/2)) {qtok2 += (nk/2);} else {qtok2 -= (nk/2);}
							if (qtok1 < (nk/2)) {qtok1 += (nk/2);} else {qtok1 -= (nk/2);}
							qtokindex[index] = (qtok3 + qtok2*nkz + qtok1*nkz*nk);




				}

			// Setup linear omega-mesh
			for (size_t i=0; i<nw; i++) omega[i] = wmin + (wmax - wmin) * float(i)/fmax(float(nw-1),float(1));
			// Now combine vectors
			for (size_t iq=0; iq<nq1*nq2*nq3; iq++) for (size_t iw=0; iw<nw; iw++) {
				size_t index = iw + iq * nw;
				indexToiq[index] = iq;
				indexToiw[index] = iw;
				QVec[index][0] = momenta(iq,0);
				QVec[index][1] = momenta(iq,1);
				QVec[index][2] = momenta(iq,2);
				QVec[index][3] = omega[iw];
				qtok[iq] = qtokindex[iq];
			}

			QVec.resize(numberOfQ/2);


		}



		//standard Write Chiq
		void writeChiqTxt(VectorSuscType& chi0Matrix) {
			std::ofstream os("chi0Full.txt");
			std::ofstream os2("chiRPA.txt");
			interaction<FieldType,psimag::Matrix,ConcurrencyType> rpa(param);
			int width(10);
			os.precision(width);
			os2.precision(width);
			os << std::fixed;
			os2 << std::fixed;
			os << "#nq1,nq2,nq3,nw: \n";
			os << "#" << nq1 << " , " << nq2 << " , " << nq3 << " , " << nw << "\n";
			std::vector<FieldType> q(3);
			for (size_t iq=0;iq<QVec.size();iq++) {
				q[0]=QVec[iq][0]; q[1]=QVec[iq][1]; q[2]=QVec[iq][2];
				os << q[0] << " " << q[1] << " " << q[2] << " " << QVec[iq][3] << " ";
				for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
					os << real(chi0Matrix[iq](l1,l2))<< " ";
				}
				for (size_t l1=0;l1<msize;l1++) for (size_t l2=0;l2<msize;l2++) {
					os << imag(chi0Matrix[iq](l1,l2))<< " ";
				}
     			ComplexType sus0(chi0Matrix[iq].calcSus());
				os << real(sus0) << " " << imag(sus0);
				os << "\n";
				SuscType chiRPA(param,conc);
     			rpa.calcRPAResult(chi0Matrix[iq],rpa.spinMatrix,chiRPA,q);
     			ComplexType susR(chiRPA.calcSus());
     			os2 << q[0] << " " << q[1] << " " << q[2] << " " << QVec[iq][3] << " ";
     			os2 << real(susR) << " "  << imag(susR) << " " << real(sus0) << " " << imag(sus0) << "\n";
			}
		}

		//optimized Write Chiq
		void writeChiqTxt(std::vector<std::vector<VectorSuscType> >& chi0Matrix,int nSite) {
			std::ofstream os("chi0Full.txt");
			std::ofstream os2("chiRPA.txt");
			interaction<FieldType,psimag::Matrix,ConcurrencyType> rpa(param);
			int nOS= param.nOrbSite;
			int width(10);
			os.precision(width);
			os2.precision(width);
			os << std::fixed;
			os2 << std::fixed;
			os << "#nq1,nq2,nq3,nw: \n";
			os << "#" << nq1 << " , " << nq2 << " , " << nq3 << " , " << nw << "\n";
			std::vector<FieldType> q(3);
			for (size_t iq=0;iq<QVec.size();iq++) {
				ComplexType sus0(0,0),susR(0,0);
				SuscType chipreRPA(param,conc,param.nSite);
				SuscType chiRPA(param,conc,param.nSite);
				ComplexMatrixType BigSpinMatrix(nOS*nOS*param.nSite,nOS*nOS*param.nSite);
				for (size_t i = 0; i < nSite; i++) {
			  	for (size_t j = i; j < nSite; j++) {
						sus0+=chi0Matrix[i][j][iq].calcSus(nOS*nOS);
						for (size_t k = 0; k < nOS; k++) {
							for (size_t l = 0; l < nOS; l++) {
								//if (conc.rank()==0) std::cerr << "mem swap" << std::endl;
								if (i==j) {
									BigSpinMatrix(k+i*nOS,l+j*nOS)= rpa.spinMatrix(k,l);
								}
								chipreRPA(k+i*nOS,l+j*nOS)= chi0Matrix[i][j][iq](k,l);
							}
						}

					}
				}
				for (size_t i = 0; i < nSite; i++) {
					for (size_t j = i+1; j < nSite; j++) {
						sus0+=conj(chi0Matrix[i][j][iq].calcSus(nOS*nOS));
						for (size_t k = 0; k < nOS; k++) {
							for (size_t l = 0; l < nOS; l++) {
								chipreRPA(l+j*nOS,k+i*nOS)= conj(chi0Matrix[i][j][iq](k,l));
							}
						}

					}
				}
			//	if (conc.rank()==0) std::cerr << "before rpa calc" << std::endl;
				rpa.calcRPAResult(chipreRPA,BigSpinMatrix,chiRPA,q);
			//	if (conc.rank()==0) std::cerr << "after rpa calc" << std::endl;
				susR = (chiRPA.calcSus(nOS*nOS*param.nSite));
			//	if (conc.rank()==0) std::cerr << "after rpa calcSus" << std::endl;
				q[0]=QVec[iq][0]; q[1]=QVec[iq][1]; q[2]=QVec[iq][2];
				os << q[0] << " " << q[1] << " " << q[2] << " " << QVec[iq][3] << " ";
				os << real(sus0) << " " << imag(sus0);
				os << "\n";
     		os2 << q[0] << " " << q[1] << " " << q[2] << " " << QVec[iq][3] << " ";
     		os2 << real(susR) << " "  << imag(susR) << " " << real(sus0) << " " << imag(sus0) << "\n";
			}
		}

		void writeChiSChiC(std::vector<ComplexMatrixType>& chiS,
						   std::vector<ComplexMatrixType>& chiN) {
			std::ofstream os1("chiSRPA.txt");
			std::ofstream os2("chiNRPA.txt");
			std::ofstream os3("chiCRPA.txt");
			int width(10);
			os1.precision(width);
			os2.precision(width);
			os3.precision(width);
			os1 << std::fixed;
			os2 << std::fixed;
			os3 << std::fixed;
			std::vector<FieldType> q(3);
			for (size_t iq=0;iq<numberOfQ;iq++) {
				q[0]=QVec[iq][0]; q[1]=QVec[iq][1]; q[2]=QVec[iq][2];
				ComplexType chiRPAS = calcChiSRPA(chiS[iq]);
				ComplexType chiRPAN = calcChiNematicRPA(chiN[iq]);
				ComplexType chiRPAC = calcChiSRPA(chiN[iq]);
     			os1 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3] << " , ";
     			os1 << real(chiRPAS) << ","  << imag(chiRPAS) <<  "\n";
     			os2 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3] << " , ";
     			os2 << real(chiRPAN) << ","  << imag(chiRPAN) <<  "\n";
     			os3 << q[0] << " , " << q[1] << " , " << q[2] << " , " << QVec[iq][3] << " , ";
     			os3 << real(chiRPAC) << ","  << imag(chiRPAC) <<  "\n";
			}
		}

		template<typename FieldType, typename ComplexMatrixType>
		void writeChi0Emery(const std::vector<std::vector<FieldType> >& qField,
						    const std::vector<ComplexMatrixType>& chi0,
						    const std::vector<ComplexMatrixType>& chi0_g,
						    const std::vector<ComplexMatrixType>& chi0_gg,
						    const std::string& filename) {
			if (conc.rank()==0) {
				std::ofstream os(filename.c_str());
				int width(10);
				os.precision(width);
				os << std::fixed;
				size_t nq(qField.size());
				for (size_t iq=0;iq<nq;iq++) {
					os << qField[iq][0] << " , " << qField[iq][1] << " , " << qField[iq][2] << " , " << qField[iq][3] << " , ";
					for (size_t i=0;i<chi0[0].n_row();i++) for (size_t j=0;j<chi0[0].n_col();j++)
						os << real(chi0[iq](i,j))<< " , " << imag(chi0[iq](i,j)) << " , ";
					for (size_t i=0;i<chi0_g[0].n_row();i++) for (size_t j=0;j<chi0_g[0].n_col();j++)
						os << real(chi0_g[iq](i,j))<< " , " << imag(chi0_g[iq](i,j)) << " , ";
					for (size_t i=0;i<chi0_gg[0].n_row();i++) for (size_t j=0;j<chi0_gg[0].n_col();j++) {
						os << real(chi0_gg[iq](i,j))<< " , " << imag(chi0_gg[iq](i,j));
						if (j<chi0_gg[0].n_col()) os << " , ";
					}
					// os << "\n";
				}
			}
		}
		template<typename FieldType, typename ComplexMatrixType>
		void readChi0Emery (std::vector<std::vector<FieldType> >& qField,
						    std::vector<ComplexMatrixType>& chi0,
						    std::vector<ComplexMatrixType>& chi0_g,
						    std::vector<ComplexMatrixType>& chi0_gg,
						    std::string& filename) {
			size_t nq(0);
			if (conc.rank()==0) {
				std::vector<FieldType> data;
				typedef std::complex<FieldType> ComplexType;
				loadVector(data,filename);
				size_t step(4+2*(3*3+19*3+19*19)); // 4 qField (qx,qy,qz,w), 3x3 for chi0, 19x3 for chi0_g, 19x19 for chi0_gg
				nq = data.size()/step;
				if (nq!=numberOfQ) {
					qField.resize(nq);
					chi0.resize(nq);
					chi0_g.resize(nq);
					chi0_gg.resize(nq);
					for (size_t iq=0;iq<nq;++iq) {
						qField[iq] = std::vector<FieldType>(3,0);
						chi0[iq] = ComplexMatrixType(3,3);
						chi0_g[iq] = ComplexMatrixType(19,3);
						chi0_gg[iq] = ComplexMatrixType(19,19);
					}
				}
				// for (size_t iq=0;iq<nq;++iq) qField[iq] = std::vector<FieldType>(3,0);
				std::cout << "number of Q points from file: " << nq << "\n";
				for (size_t iq=0;iq<nq;iq++) {
					qField[iq][0] = data[0 + iq*step];
					qField[iq][1] = data[1 + iq*step];
					qField[iq][2] = data[2 + iq*step];
					qField[iq][3] = data[3 + iq*step];
					size_t ishift(4 + iq*step);
					for (size_t i=0;i<chi0[iq].n_row();i++) for (size_t j=0;j<chi0[iq].n_col();j++) {
						FieldType r1 = data[ishift + 2*(j+i*chi0[0].n_col())];
						FieldType i1 = data[ishift + 2*(j+i*chi0[0].n_col()) + 1];
						chi0[iq](i,j) = ComplexType(r1,i1);
					}
					ishift = 4 + iq*step + 2*chi0[0].n_row()*chi0[0].n_col();
					for (size_t i=0;i<chi0_g[iq].n_row();i++) for (size_t j=0;j<chi0_g[iq].n_col();j++) {
						FieldType r1 = data[ishift + 2*(j+i*chi0_g[0].n_col())];
						FieldType i1 = data[ishift + 2*(j+i*chi0_g[0].n_col()) + 1];
						chi0_g[iq](i,j) = ComplexType(r1,i1);
					}
					ishift = 4 + iq*step + 2*chi0[0].n_row()*chi0[0].n_col() + 2*chi0_g[0].n_row()*chi0_g[0].n_col();
					for (size_t i=0;i<chi0_gg[iq].n_row();i++) for (size_t j=0;j<chi0_gg[iq].n_col();j++) {
						FieldType r1 = data[ishift + 2*(j+i*chi0_gg[0].n_col())];
						FieldType i1 = data[ishift + 2*(j+i*chi0_gg[0].n_col()) + 1];
						chi0_gg[iq](i,j) = ComplexType(r1,i1);
					}
				}
			} else {
				conc.broadcast(nq);
				if (nq!=numberOfQ) {
					qField.resize(nq);
					chi0.resize(nq);
					chi0_g.resize(nq);
					chi0_gg.resize(nq);
					for (size_t iq=0;iq<nq;++iq) {
						qField[iq] = std::vector<FieldType>(3,0);
						chi0[iq] = ComplexMatrixType(3,3);
						chi0_g[iq] = ComplexMatrixType(19,3);
						chi0_gg[iq] = ComplexMatrixType(19,19);
					}
				}
			}
			numberOfQ = nq;
			for (size_t iq=0;iq<numberOfQ;iq++) {
				conc.broadcast(qField[iq]);
				conc.broadcast(chi0[iq]);
				conc.broadcast(chi0_g[iq]);
				conc.broadcast(chi0_gg[iq]);
			}
		}


		void printGap(const std::string& file="gapRPA.txt") {
			VectorType data;
			loadVector(data,file);
			size_t step=5;
			std::cout << "File=" << file << "\n";
			size_t nk(data.size()/step);
			std::cout << "number of k-points in gap file = " << nk << "\n";

			GapType Delta(param,conc);
			// momentumDomain<Field,psimag::Matrix> kmesh(param,param.nkInt,param.nkIntz,param.dimension);
			// kmesh.set_momenta(false);
			// BandsType bands(param,conc,kmesh,true);
			// VectorType ek(param.nOrb,0), k(3);
			// ComplexMatrixType ak(param.nOrb,param.nOrb);
			std::ofstream os("gap.txt");
			int width(10);
			os.precision(width);
			os << std::fixed;
			for (size_t ik = 0; ik < nk; ++ik)	{
				// kmesh.momenta.getRow(ik,k);
				// bands.getEkAndAk(k,ek,ak);
				VectorType k(3,0); size_t iband;
				k[0] = data[ik*step];
				k[1] = data[ik*step+1];
				k[2] = data[ik*step+2];
				iband = size_t(data[ik*step+3]);
				FieldType gapRPA(data[ik*step+4]);
				// for (size_t iband=0;iband<param.nOrb;iband++) {
				ComplexType gap1 = Delta(k,iband);
					// gap1 *= pow(param.Omega0,2)/(pow(ek[iband],2)+pow(param.Omega0,2)); // Lorentzian cut-off

				os << k[0] << " , " << k[1] << " , " << k[2] << " , "
					   << iband << " , " << gapRPA << " , " << real(gap1) << "\n";
				std::cout << k[0] << " , " << k[1] << " , " << k[2] << " , "
					   << iband << " , " << gapRPA << " , " << real(gap1) << "\n";
				// }

			}
		}

		void printGap2() {
			GapType Delta(param,conc);
			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			kmesh.set_momenta(false);
			BandsType bands(param,conc,kmesh,true);
			VectorType ek(param.nOrb,0), k(3);
			ComplexMatrixType ak(param.nOrb,param.nOrb);
			std::ofstream os("gapAll.txt");
			int width(10);
			os.precision(width);
			os << std::fixed;
			for (size_t ik = 0; ik < kmesh.nktot; ++ik)	{
				kmesh.momenta.getRow(ik,k);
				bands.getEkAndAk(k,ek,ak);
				std::vector<FieldType> gap1(param.nOrb);
				for (size_t iband=0;iband<param.nOrb;iband++) {
					gap1[iband] = real(Delta(k,iband));
					gap1[iband] *= pow(param.Omega0,2)/(pow(ek[iband],2)+pow(param.Omega0,2)); // Lorentzian cut-off
				}

				os << k[0] << "  " << k[1] << "  " << k[2] << "  "
					   << gap1 << "\n";
			}
		}

		void printGap3() {
			ferminator<FieldType,BandsType,psimag::Matrix,ConcurrencyType> FSpoints(param,conc,1);
			size_t nk(FSpoints.nTotal);
			momentumDomain<Field,psimag::Matrix,ConcurrencyType> kmesh(param,conc,param.nkInt,param.nkIntz,param.dimension);
			BandsType bands(param,conc,kmesh,true);
			VectorType ek(param.nOrb,0);
			ComplexMatrixType ak(param.nOrb,param.nOrb);
			GapType Delta(param,conc);
			std::ofstream os("gapOnFS.txt");
			int width(10);
			os.precision(width);
			os << std::fixed;
			for (size_t ik = 0; ik < nk; ++ik)	{
				VectorType k(3,0); size_t iband;
				k[0] = FSpoints.kFx[ik];
				k[1] = FSpoints.kFy[ik];
				k[2] = FSpoints.kFz[ik];
				iband = FSpoints.kFtoBand[ik];
				bands.getEkAndAk(k,ek,ak);
				ComplexType gap1 = Delta(k,iband,ak);
				// gap1 *= pow(param.Omega0,2)/(pow(ek[iband],2)+pow(param.Omega0,2)); // Lorentzian cut-off
				os << k[0] << " , " << k[1] << " , " << k[2] << " , "
					   << iband;
				for (size_t iorb=0;iorb<param.nOrb;iorb++) os << " , " << FSpoints.weights[ik][iorb];
			    os << " , " << real(gap1) << "\n";
			}
		}


	};
}

#endif
