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
    double Max, Min, diffMax, err, maxpw_relerr = 0, relerr, normErr = 0;
    Max = oriData[0];
    Min = oriData[0];
    diffMax = fabs(deData[0] - oriData[0]);
	
    for (i = 0; i < length; i++)
    {
    	if (Max < oriData[i]) Max = oriData[i];
    	if (Min > oriData[i]) Min = oriData[i];
		err = fabs(deData[i] - oriData[i]);
		normErr += err*err;
    	if (diffMax < err)
    		diffMax = err;
        if(oriData[i]!=0)
        {
            relerr = err/fabs(oriData[i]);
            if(maxpw_relerr<relerr)
                maxpw_relerr = relerr;
        }
    }
	std::cout << "check: " << normErr << std::endl;
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
	std::cout << "Number of cells in the file: " << num_cells << std::endl;

	double* xx = new double[num_cells];
	double* yy = new double[num_cells];
	double* zz = new double[num_cells];
    double* oriVol = new double[num_cells];
	double* vol = new double[num_cells];
	double* rho = new double[num_cells];
	double* temp = new double[num_cells];

	size_t rate_2 = 0;

    double maxVol = -DBL_MAX;
    double minVol = DBL_MAX;
	double maxRho = -DBL_MAX;
    double minRho = DBL_MAX;

	for (size_t i = 0; i < num_cells; ++i) 
	{
		f.read( reinterpret_cast<char*>(&xx[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&yy[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&zz[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&oriVol[i]) , sizeof(double));
		f.read( reinterpret_cast<char*>(&rho[i]) , sizeof(double));
        if (rho[i] > maxRho)
			maxRho = rho[i];
		if (rho[i] < minRho)
			minRho = rho[i];
		f.read( reinterpret_cast<char*>(&temp[i]) , sizeof(double));
        vol[i] = std::cbrt(oriVol[i]);
        if (vol[i] > maxVol)
			maxVol = vol[i];
		if (vol[i] < minVol)
			minVol = vol[i];
	}

	double range = maxRho - minRho;
	std::cout << "range:" << range << std::endl;

    size_t grid_0 = 512;
	size_t grid_2 = 256;

    char orderName[640];
	sprintf(orderName, "%s", argv[2]);
    std::ifstream f_2;
    f_2.open(orderName, std::ifstream::binary);
	if (f_2.fail()) {
		std::cout << "Error opening file" << std::endl;
		return 0;
	}

    double* reorder = new double[grid_2 * grid_2];
    double* degrid_0 = new double[grid_0 * grid_0 * grid_0];
    double* degrid_2 = new double[grid_2 * grid_2 * grid_2];
    double* densgrid_0 = new double[grid_0 * grid_0 * grid_0];
    double* densgrid_2 = new double[grid_2 * grid_2 * grid_2];
    int* mark_0 = new int[grid_0 * grid_0 * grid_0];
    int* mark_2 = new int[grid_2 * grid_2 * grid_2];
    std::fill(mark_0, mark_0+grid_0 * grid_0 * grid_0, -1);
    std::fill(mark_2, mark_2+grid_2 * grid_2 * grid_2, -1);
    
    double chechSum_all = 0;

    for (size_t i = 0; i < num_cells; ++i)
	{
		size_t xi, yi, zi, cpos, mx, my, mz;
		// Debug: Make sure 'round' is working right
		double dxi = (xx[i] / minVol);
		if (round(dxi) != dxi)
		{
			std::cout.precision(10);
			if( abs(dxi - round(dxi)) > 0.5)
				cout << "WARNING: Potential mapping error: " << dxi << " vs " << round(dxi) << std::endl;
			
		}

		xi = round(xx[i] / minVol);
		yi = round(yy[i] / minVol);
		zi = round(zz[i] / minVol);

		if (vol[i] == minVol)
		{
			cpos = xi + (yi * (grid_0)) + (zi * (grid_0) * (grid_0));
			densgrid_0[cpos] = rho[i];
		}
		else
		{
			rate_2 = vol[i]/minVol;	
            cpos = xi/rate_2 + (yi/rate_2 * (grid_2)) + (zi/rate_2 * (grid_2) * (grid_2));
            densgrid_2[cpos] = rho[i];
		}
	}

    for (size_t i = 0; i < grid_2 * grid_2; ++i) 
	{
        f_2.read( reinterpret_cast<char*>(&reorder[i]) , sizeof(double));
    }


    double* coarse_z = new double[num_cells];

    size_t cnt = 0;
    for (size_t z = 0; z < grid_2; ++z) {
        for (size_t y = 0; y < grid_2; ++y) {
            for (size_t x = 0; x < grid_2; ++x) {
                if (densgrid_2[int(reorder[x + y * grid_2]) + z * grid_2 * grid_2] != 0 ){

                    coarse_z[x + y * grid_2 + z * grid_2* grid_2 + cnt] = densgrid_2[int(reorder[x + y * grid_2]) + z * grid_2 * grid_2];
                    mark_2[int(reorder[x + y * grid_2]) + z * grid_2 * grid_2] = x + y * grid_2 + z * grid_2* grid_2 + cnt;
                } else {
                    size_t j = int(floor(reorder[x + y * 256] / 256));
                    size_t i = int(reorder[x + y * 256] - j * 256);
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt] = densgrid_0[i * 2 + (j * 2 * 512) + (z * 2 * grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 1] = densgrid_0[(i * 2 + 1) + (j * 2 * 512) + (z * 2 * grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 2] = densgrid_0[(i * 2) + ((j * 2 + 1) * 512) + (z * 2 * grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 3] = densgrid_0[(i * 2 + 1) + ((j * 2 + 1) * 512) + (z * 2 * grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 4] = densgrid_0[i * 2 + (j * 2 * 512) + ((z * 2 + 1)* grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 5] = densgrid_0[(i * 2 + 1) + (j * 2 * 512) + ((z * 2 + 1)* grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 6] = densgrid_0[(i * 2) + ((j * 2 + 1) * 512) + ((z * 2 + 1)* grid_0 * grid_0)];
                    coarse_z[x + y * 256 + z * grid_2 * grid_2 + cnt + 7] = densgrid_0[(i * 2 + 1) + ((j * 2 + 1) * 512) + ((z * 2 + 1)* grid_0 * grid_0)];

                    mark_0[i * 2 + (j * 2 * 512) + (z * 2 * grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt;
                    mark_0[(i * 2 + 1) + (j * 2 * 512) + (z * 2 * grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt + 1;
                    mark_0[(i * 2) + ((j * 2 + 1) * 512) + (z * 2 * grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt + 2;
                    mark_0[(i * 2 + 1) + ((j * 2 + 1) * 512) + (z * 2 * grid_0 * grid_0)]= x + y * 256 + z * grid_2 * grid_2 + cnt + 3;
                    mark_0[i * 2 + (j * 2 * 512) + ((z * 2 + 1)* grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt + 4;
                    mark_0[(i * 2 + 1) + (j * 2 * 512) + ((z * 2 + 1)* grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt + 5;
                    mark_0[(i * 2) + ((j * 2 + 1) * 512) + ((z * 2 + 1)* grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt + 6;
                    mark_0[(i * 2 + 1) + ((j * 2 + 1) * 512) + ((z * 2 + 1)* grid_0 * grid_0)] = x + y * 256 + z * grid_2 * grid_2 + cnt + 7;

                    cnt += 7;
                }
            }
        }
    }
    

    char *cfg = argv[3];
    double eb = atof(argv[4]);
    char path[650];
    sprintf(path, "%s.zmesh", inName);
    size_t outSize = 0;
    coarse_z = AMR_compress(coarse_z, cfg, num_cells, eb, path, &outSize);
    std::cout << "zMesh compression ratio: " << double(num_cells) * 8 / outSize << std::endl;
	std::cout << "zMesh bitrate: " << double(outSize) * 8 / num_cells << std::endl;

    for (size_t i = 0; i < grid_2 * grid_2 * grid_2; ++i) 
	{
        if (mark_2[i] != -1) {
            degrid_2[i] = coarse_z[mark_2[i]];
        }
    }
    for (size_t i = 0; i < grid_0 * grid_0 * grid_0; ++i) 
	{
        if (mark_0[i] != -1) {
            degrid_0[i] = coarse_z[mark_0[i]];
        }
    }

    double normErr_0 = 0, normErr_2 = 0, mse = 0, mse_0 = 0 , mse_2 = 0;
	for (size_t i = 0; i < grid_0*grid_0*grid_0; ++i) 
	{
		if (densgrid_0[i] != 0)
		{
			// std::cout << origrid_0[i] << " " << densgrid_0[i] << std::endl;
			normErr_0 += (densgrid_0[i] - degrid_0[i]) * (densgrid_0[i] - degrid_0[i]) ;
		}
	}
	for (size_t i = 0; i < grid_2*grid_2*grid_2; ++i) 
	{
		if (densgrid_2[i] != 0)
		{
			normErr_2 += (densgrid_2[i] - degrid_2[i]) * (densgrid_2[i] - degrid_2[i]);
		}
	}

    mse = (normErr_0  + pow(2,3)*normErr_2)/(grid_0 * grid_0 * grid_0);
    std::cout << "zMesh PSNR: " << 20*log10(range)-10*log10(mse) << std::endl;

    return 0;
}

