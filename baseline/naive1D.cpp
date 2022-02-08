using namespace std;
#include <cfloat>
#include <cmath>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>		// std::string
#include <sstream>		// std::stringstream
#include <vector>
#include <map>
#include <algorithm>
#include "sz.h"
#include "rw.h"

double *AMR_compress(double *oriData, char* cfgFile, size_t length, double eb, char* zipFilePath, size_t* comSize)
{	
    int status = 0;
    printf("cfgFile=%s\n", cfgFile); 
    status = SZ_Init(cfgFile);
	
    size_t outSize = 0;
	// std::cout << length4d << std::endl;
	// std::cout << blksize << std::endl;
	unsigned char* bytes = SZ_compress_args(SZ_DOUBLE, oriData, &outSize, ABS, eb, 1, 1, 0, 0, 0, 0, length);
    // unsigned char *bytes = SZ_compress(SZ_DOUBLE, oriData, &outSize, 0, length4d, blksize, blksize, blksize);
    writeByteData(bytes, outSize, zipFilePath, &status);
	if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be written!\n", zipFilePath);
	exit(0);
    }
	std::cout << "compressedFile=" << zipFilePath << std::endl;
	*comSize += outSize;
	double *deData = (double *)SZ_decompress(SZ_DOUBLE, bytes, outSize, 0, 0, 0, 0, length);

 	size_t i;
    double Max, Min, diffMax, err, maxpw_relerr = 0, relerr;
    Max = oriData[0];
    Min = oriData[0];
    diffMax = fabs(deData[0] - oriData[0]);
	
    for (i = 0; i < length; i++)
    {
    	if (Max < oriData[i]) Max = oriData[i];
    	if (Min > oriData[i]) Min = oriData[i];
		err = fabs(deData[i] - oriData[i]);
    	if (diffMax < err)
    		diffMax = err;
        if(oriData[i]!=0)
        {
            relerr = err/fabs(oriData[i]);
            if(maxpw_relerr<relerr)
                maxpw_relerr = relerr;
        }
    }
    printf ("Max absolute error = %.20G\n", diffMax);
    printf ("Max relative error = %.20G\n", diffMax/(Max-Min));
    printf ("Max pw_relative err = %.20G\n", maxpw_relerr);

    if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be written!\n");
	// free(oriData);
	exit(0);
    }
	// delete[] oriData;
	delete[] bytes;
    // free(oriData);
    // free(bytes);
    SZ_Finalize();
	printf("done\n");
	return deData;
    
}

int main(int argc, char *argv[]) {
    char inName[640];
	sprintf(inName, "%s", argv[1]);
    std::ifstream f;
    f.open(inName, std::ifstream::binary);
	if (f.fail()) {
		std::cout << "Error opening file" << std::endl;
		return 0;
	}

    f.seekg(0, std::ios::end);
	size_t file_len = f.tellg();
	f.seekg(0, std::ios::beg);

	size_t num_cells = file_len / (6 * sizeof(double)); // changes depending on # of fields
	// num_cells = 60000000;
	std::cout << "Number of cells in the file: " << num_cells << std::endl;

	double* xx = new double[num_cells];
	double* yy = new double[num_cells];
	double* zz = new double[num_cells];
	double* vol = new double[num_cells];
	double* oriVol = new double[num_cells];
	double* rho = new double[num_cells];
	double* temp = new double[num_cells];
    size_t size[4] = {0, 0, 0, 0};
	size_t rate[4] = {1, 0, 0, 0};
	double maxRho = -DBL_MAX;
    double minRho = DBL_MAX;

	double lVoltemp[4];
	size_t iden = 0;
	bool add = 1;

	for (size_t i = 0; i < num_cells; ++i) 
	{
		f.read( reinterpret_cast<char*>(&xx[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&yy[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&zz[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&oriVol[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&rho[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&temp[i]) , sizeof(double));
		vol[i] = std::cbrt(oriVol[i]);
		if (rho[i] > maxRho)
			maxRho = rho[i];
		if (rho[i] < minRho)
			minRho = rho[i];

        if (iden == 0)
		{
			lVoltemp[0] = vol[i];
			iden += 1;
		} else {
			for (size_t j = 0; j < iden; j++){
				if (vol[i] == lVoltemp[j]){
					add = 0;
				}			
			}	

			if (add == 1){
				// std::cout << iden << std::endl;
				lVoltemp[iden] = vol[i];
				iden += 1;
			}
			add = 1;
		}
	}

	double range = maxRho - minRho;
	std::cout << "range : " << range  << std::endl;
	double lVol[iden];
	std::cout << "number of AMR levels: " << iden << std::endl;

	double **temprho = new double*[iden];
    for(int i = 0; i < iden; ++i) {
        temprho[i] = new double[num_cells]{0.0};
    }	

	for (size_t i = 0; i < iden; i++)
	{
		lVol[i] = lVoltemp[i];
	}
	std::sort(lVol, lVol + iden);

	for (size_t i = 0; i < num_cells; ++i){
        if (vol[i] == lVol[0]){
            temprho[0][size[0]] = rho[i];
            size[0]++;
        }
        else if (vol[i] == lVol[1]){
			rate[1] = vol[i]/lVol[0];
            temprho[1][size[1]] = rho[i];
            size[1]++;
        }    
		else if (vol[i] == lVol[2])
		{
			rate[2] = vol[i]/lVol[0];
            temprho[2][size[2]] = rho[i];
            size[2]++;
		}
        else if (vol[i] == lVol[3]){
			rate[3] = vol[i]/lVol[0];
            temprho[3][size[3]] = rho[i];
            size[3]++;
        }        
    }
    
	double **lRho = new double*[iden];
    for(size_t i = 0; i < iden; ++i) {
        lRho[i] = new double[size[i]]{0.0};
    }	
	for(size_t j = 0; j < iden; ++j) {
        for(size_t i = 0; i < size[j]; ++i) {
			lRho[j][i] = temprho[j][i];
		}
    }

	double **de = new double*[iden];
    for(size_t i = 0; i < iden; ++i) {
        de[i] = new double[size[i]]{0.0};
    }	

    char *cfg = argv[2];
	double eb[4] = {atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6])};

	size_t outSize = 0;

	for(size_t i = 0; i < iden; ++i) {
		char path[650];
		sprintf(path, "%s.1d_%d", inName, i);
        de[i] = AMR_compress(lRho[i], cfg, size[i], eb[i], path, &outSize);
    }	

	std::cout << "naive1D compression ratio: " << double(num_cells) * 8 / outSize << std::endl;
	std::cout << "naive1D bitrate: " << double(outSize) * 8 / num_cells << std::endl;

	double normErr[iden] = {0, 0, 0, 0};

	for (size_t j = 0; j < iden; ++j) {
		for (size_t i = 0; i < size[j]; ++i) 
		{
			normErr[j] += (de[j][i] - lRho[j][i]) * (de[j][i] - lRho[j][i]);
		}
	}

	double mse = 0;
	for (size_t i = 0; i < iden; ++i) {
		// std::cout << rate[i] << std::endl;
		// std::cout << normErr[i] << std::endl;
		mse += pow(rate[i],3) * normErr[i];
	}

	size_t sAll = 0;
	for (size_t i = 0; i < iden; ++i) {
		sAll += pow(rate[i],3) * size[i];
	}
	mse = mse / sAll;

	std::cout << "naive1D PSNR: " << 20*log10(range)-10*log10(mse) << std::endl;

    return 0;
}

