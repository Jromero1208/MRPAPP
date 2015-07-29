//-*-C++-*-

#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "parameters.h"
#include "momentumDomain.h"
#include "bandstructure.h"

namespace rpa {

	template<typename Field, template<typename> class MatrixTemplate, typename ConcurrencyType>
	class greensFunction:
		public MatrixTemplate<std::complex<Field> >
 	{

	private:
		typedef Field 				                     FieldType;
		typedef std::complex<Field> 							 ComplexType;
		typedef psimag::Matrix<ComplexType> 			 ComplexMatrixType;
		typedef std::vector<Field> 								 VectorType;

		const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& param;
		size_t nOrb;
		size_t nktot;
		ConcurrencyType& concurrency;
		VectorType k1,k2,k3,w;

 	public:
 		typedef psimag::Matrix<std::complex<Field> > BaseType;
		momentumDomain<Field,psimag::Matrix,ConcurrencyType> momentumDomain1;
		size_t nwn;
 		BaseType& super;


	greensFunction(const rpa::parameters<Field,MatrixTemplate,ConcurrencyType>& parameters,
		ConcurrencyType& conc):
		BaseType(size_t(pow(parameters.nkInt,parameters.dimension))*parameters.nwn,parameters.nOrb*parameters.nOrb),
		param(parameters),
		nOrb(param.nOrb),
		nktot(size_t(pow(param.nkInt,param.dimension))),
		concurrency(conc),
		momentumDomain1(param,concurrency,param.nkInt,param.dimension),
		nwn(param.nwn),
		super (*this)
		{

			if (!param.Green) return;
			if (param.Gfile != "")
			{
				readCSVFile();
			} else {
			momentumDomain1.set_momenta(true);
			}
		}


	void toZero() {
		size_t dim1(nktot*param.nwn);
		size_t dim2(nOrb*nOrb);
		for (size_t i=0; i<dim1; i++) {
			for (size_t j=0; j<dim2; j++) {
				(*this)(i,j) = ComplexType(0.0);
			}
		}
	}


	void init() {

		if (param.Gfile != "") return;

		std::cout << "Starting initialization of G(k,iw) ... \n";
		size_t nktot(momentumDomain1.nktot);
		VectorType ek(nOrb);
		ComplexMatrixType ak(nOrb,nOrb);
		VectorType k(3),q(3),kq(3);
		bandstructure<Field,psimag::Matrix,ConcurrencyType> bandstructure(param,concurrency);


		for (size_t ik = 0; ik < nktot; ++ik)
		{
				momentumDomain1.momenta.getRow(ik,k);

				bandstructure.getEkAndAk(k,ek,ak);

				for (int in = 0; in < nwn; ++in)
				{
					if (ik == 0) w.push_back((2*in+1)*param.pi_f);
					Field iwn = (2*in+1)*param.pi_f*param.temperature;
					size_t ikw = in + (ik * nwn);
					for (size_t l1 = 0; l1 < nOrb; ++l1)
					{
						for (size_t l2 = 0; l2 < nOrb; ++l2)
						{
							size_t ind = l2 + l1 * nOrb;
							ComplexType c2(0.0);
							for (size_t iband = 0; iband < nOrb; ++iband)
							{
								ComplexType c1 = ak(l1,iband) * conj(ak(l2,iband));
								c2 += ComplexType(1.)/(ComplexType(-ek[iband],iwn)) * c1;
							}
							(*this)(ikw,ind) = c2;
						}
					}
				}
			//		std::cout << k << " " << ek << " " << ak << std::endl;
		}
		std::cout << "G(k,iw) initialized... \n";
	}




	void readCSVFile() {
		std::string file = param.Gfile;

		std::cout << "Reading Green's Function from file: " << file << std::endl;
		VectorType data;
		loadVector(data,file);
		// We assume that each line has the format k1 k2 k3 w (ReG00,ImG00) (ReG01,ImG01)... (ReGMM,ImGMM)
		size_t length = 4 + (2 * param.nOrb * param.nOrb);
		size_t nLinesTotal(data.size()/length);
		if (concurrency.rank()==0) std::cout << "Green's Function file contains " << nLinesTotal << " lines\n";
		for (size_t i = 0; i < nLinesTotal; i++) {

			if (data[i*length+3] == data[3]) {
			  k1.push_back(data[i*length]);
				k2.push_back(data[i*length+1]);
				k3.push_back(data[i*length+2]);
				//std::cerr << "i, k1 " << i << ", " << k1[i*length] << std::endl;
			}

			if (i == 0 ||  w[w.back()] < data[i*length+3])  w.push_back(data[i*length+3]);
			for (size_t j = 0; j < (param.nOrb * param.nOrb); j++) {
				(*this)(i,j) = ComplexType(data[i*length+4+2*j],data[i*length+5+2*j]);
			}
		}
		data.resize(0);

		momentumDomain1.set_momenta(k1,k2,k3);

		size_t nLines = k1.size();
		if (concurrency.rank()==0) std::cout << nLines <<" Total Kpoints used from Green's Function input file\n";
	}

	void printGreens() {

		std::ofstream os("GreensFunction.txt");
		int width(10);
		os.precision(width);
		os << std::fixed;
//		os << "# nk1, nk2, nk3, nw: \n";
//		os << "# "<< nk2 << " , " << nk2 << " , " << nk3 << " , " << nw << "\n";
		std::vector<FieldType> k(3);
		for (size_t i=0; i < nktot;i++) {
			k[0]=momentumDomain1.momenta[i][0]; k[1]=momentumDomain1.momenta[i][1]; k[2]=momentumDomain1.momenta[i][2];
			for (size_t j=0; j < nwn; j++) {
				os << k[0] << " , " << k[1] << " , " << k[2] << " , " << w[j];
				for (size_t k = 0; k < (param.nOrb * param.nOrb); k++) {
					os << " , " << real((*this)(j+i*nwn,k)) << " , " << imag((*this)(j+i*nwn,k));
				}
			os << "\n";
			}
		}
	}

	};
}

#endif
