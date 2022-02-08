using namespace std;
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>		// std::string
#include <sstream>		// std::stringstream
#include <vector>
#include <map>
#include <algorithm>
#include "sz.h"
#include "rw.h"
#include <chrono>
#include <bits/stdc++.h>

double *AMR_compress(double *oriData, char* cfgFile, size_t blksize_x, size_t blksize_y, size_t blksize_z, size_t length4d, double eb, char* zipFilePath, size_t* comSize)
{	
    int status = 0;
    // printf("cfgFile=%s\n", cfgFile); 
    status = SZ_Init(cfgFile);
	
    size_t outSize = 0;
	// std::cout << length4d << std::endl;
	// std::cout << blksize << std::endl;
	unsigned char* bytes = SZ_compress_args(SZ_DOUBLE, oriData, &outSize, ABS, eb, 1, 1, 0, length4d, blksize_z, blksize_y, blksize_x);
    // unsigned char *bytes = SZ_compress(SZ_DOUBLE, oriData, &outSize, 0, length4d, blksize, blksize, blksize);
    writeByteData(bytes, outSize, zipFilePath, &status);
	if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be written!\n", zipFilePath);
	exit(0);
    }
	std::cout << "compressedFile=" << zipFilePath << std::endl;
	*comSize += outSize;
	double *deData = (double *)SZ_decompress(SZ_DOUBLE, bytes, outSize, 0, length4d, blksize_z, blksize_y, blksize_x);

 	size_t i;
    double Max, Min, diffMax, err, maxpw_relerr = 0, relerr;
    Max = oriData[0];
    Min = oriData[0];
    diffMax = fabs(deData[0] - oriData[0]);

	if (length4d == 0){
		length4d = 1;
	}
	
    for (i = 0; i < length4d*blksize_x*blksize_y*blksize_z; i++)
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

double *Opst_compress(double *densgrid, char* tp3d, size_t grid, size_t blkSize, char* inName, double eb, char* cfg, size_t lvl, size_t* comSize)
{

    size_t blkNum = grid / blkSize;
    size_t* big = new size_t[blkNum*blkNum*blkNum]{0};
    size_t maxBlk =0;
    for (size_t z = 0; z < blkNum; z++){
		for (size_t y = 0; y < blkNum; y++){
			for (size_t x = 0; x < blkNum; x++){
                if (tp3d[x + (y * blkNum) + (z * blkNum * blkNum)] == 1) {
                    if (x == 0 || y == 0 || z == 0) {
                        big[x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
                    }else {
                        big[x + (y * blkNum) + (z * blkNum * blkNum)] = min(min(min(min(min(min(big[(x-1) + ((y-1) * blkNum) + (z * blkNum * blkNum)], big[(x - 1) + (y * blkNum) + (z * blkNum * blkNum)]), big[x + ((y-1) * blkNum) + (z * blkNum * blkNum)]), big[x + (y * blkNum) + ((z-1) * blkNum * blkNum)]), big[(x-1) + (y * blkNum) + ((z-1) * blkNum * blkNum)]), big[x + ((y-1) * blkNum) + ((z-1) * blkNum * blkNum)]), big[(x-1) + ((y-1) * blkNum) + ((z-1) * blkNum * blkNum)]) + 1;
                    }
                }
                maxBlk = max(maxBlk, big[x + (y * blkNum) + (z * blkNum * blkNum)]);
            }
        }
    }
    // std::cout << "max: " << maxBlk << std::endl;

    double **ary = new double*[maxBlk];
    for(int i = 0; i < maxBlk; ++i) {
        ary[i] = new double[grid*grid*grid];
    }

	char * take = new char[blkNum*blkNum*blkNum]{0};
	size_t sumBlk = 0;
    size_t * c = new size_t[maxBlk]{0};
    size_t * ind = new size_t[maxBlk]{0};
	for (int z = blkNum - 1; z > -1; z--){
		for (int y = blkNum - 1; y > -1; y--){
			for (int x = blkNum - 1; x > -1; x--){
				size_t size = big[x + (y * blkNum) + (z * blkNum * blkNum)];
				take[x + (y * blkNum) + (z * blkNum * blkNum)] = size;
				if (size > 0)
				{
					sumBlk += size * size * size;
					for (size_t i = z - size + 1; i < z + 1; i++){
						for (size_t j = y - size + 1; j < y + 1; j++){
							for (size_t k = x - size + 1; k < x + 1; k++){
								big[k + (j * blkNum) + (i * blkNum * blkNum)] = 0;
								tp3d[k + (j * blkNum) + (i * blkNum * blkNum)] -= 1;
							}
						}
					}

                    for (size_t k = (z - size + 1) * blkSize; k < (z + 1) * blkSize; k++){
                        for (size_t j = (y - size + 1) * blkSize; j < (y + 1) * blkSize; j++){
                            for (size_t i = (x - size + 1) * blkSize; i < (x + 1)  * blkSize; i++){
                                ary[size - 1][ind[size - 1]] = densgrid[i + (j * grid) + (k * grid * grid)];
                                ind[size - 1] ++;
                            }
                        }
                    }
                    c[size - 1]++;

					size_t xb = min(x+maxBlk, blkNum);
					size_t yb = min(y+maxBlk, blkNum);
					for (size_t k = z - size + 1; k < z + 1; k++){
						for (size_t j = y - size + 1; j < y + 1; j++){
							for (size_t i = x + 1; i < xb; i++){
								if (tp3d[i + (j * blkNum) + (k * blkNum * blkNum)] == 1) {
									if (i == 0 || j == 0 || k == 0) {
										big[i + (j * blkNum) + (k * blkNum * blkNum)] = 1;
									}else {
										big[i + (j * blkNum) + (k * blkNum * blkNum)] = min(min(min(min(min(min(big[(i-1) + ((j-1) * blkNum) + (k * blkNum * blkNum)], big[(i - 1) + (j * blkNum) + (k * blkNum * blkNum)]), big[i + ((j-1) * blkNum) + (k * blkNum * blkNum)]), big[i + (j * blkNum) + ((k-1) * blkNum * blkNum)]), big[(i-1) + (j * blkNum) + ((k-1) * blkNum * blkNum)]), big[i + ((j-1) * blkNum) + ((k-1) * blkNum * blkNum)]), big[(i-1) + ((j-1) * blkNum) + ((k-1) * blkNum * blkNum)]) + 1;
									}
								}
							}
						}
					}

					for (size_t k = z - size + 1; k < z + 1; k++){
						for (size_t j = y + 1; j <  yb; j++){
							for (size_t i = x - size + 1; i < x + 1; i++){
								if (tp3d[i + (j * blkNum) + (k * blkNum * blkNum)] == 1) {
									if (i == 0 || j == 0 || k == 0) {
										big[i + (j * blkNum) + (k * blkNum * blkNum)] = 1;
									}else {
										big[i + (j * blkNum) + (k * blkNum * blkNum)] = min(min(min(min(min(min(big[(i-1) + ((j-1) * blkNum) + (k * blkNum * blkNum)], big[(i - 1) + (j * blkNum) + (k * blkNum * blkNum)]), big[i + ((j-1) * blkNum) + (k * blkNum * blkNum)]), big[i + (j * blkNum) + ((k-1) * blkNum * blkNum)]), big[(i-1) + (j * blkNum) + ((k-1) * blkNum * blkNum)]), big[i + ((j-1) * blkNum) + ((k-1) * blkNum * blkNum)]), big[(i-1) + ((j-1) * blkNum) + ((k-1) * blkNum * blkNum)]) + 1;
									}
								}
							}
						}
					}

					for (size_t k = z - size + 1; k < z + 1; k++){
						for (size_t j = y + 1; j < yb; j++){
							for (size_t i = x + 1; i < xb; i++){
								if (tp3d[i + (j * blkNum) + (k * blkNum * blkNum)] == 1) {
									if (i == 0 || j == 0 || k == 0) {
										big[i + (j * blkNum) + (k * blkNum * blkNum)] = 1;
									}else {
										big[i + (j * blkNum) + (k * blkNum * blkNum)] = min(min(min(min(min(min(big[(i-1) + ((j-1) * blkNum) + (k * blkNum * blkNum)], big[(i - 1) + (j * blkNum) + (k * blkNum * blkNum)]), big[i + ((j-1) * blkNum) + (k * blkNum * blkNum)]), big[i + (j * blkNum) + ((k-1) * blkNum * blkNum)]), big[(i-1) + (j * blkNum) + ((k-1) * blkNum * blkNum)]), big[i + ((j-1) * blkNum) + ((k-1) * blkNum * blkNum)]), big[(i-1) + ((j-1) * blkNum) + ((k-1) * blkNum * blkNum)]) + 1;
									}
								}
							}
						}
					}
				}
            }
        }
    }


    // for(int i = 0; i < maxBlk; ++i) {
    //     std::cout << c[i] << " ";
    // }
    // std::cout << std::endl;
    // for(int i = 0; i < maxBlk; ++i) {
    //     std::cout << ind[i] << " ";
    // }
    // std::cout << std::endl;

    char path[maxBlk][650];

    double **ts = new double*[maxBlk];
    for(int i = 0; i < maxBlk; ++i) {
        sprintf(path[i], "%s.3d_%d_%d", inName,lvl , i+1);
        ts[i] = new double[ind[i]];
        for (size_t j = 0; j < ind[i]; j++)
        {
            ts[i][j] = ary[i][j];
        }
        if (ind[i] > 0) {
			if (c[i] == 1){
				c[i] = 0;
			}
            ts[i] = AMR_compress(ts[i], cfg, blkSize * (i + 1), blkSize * (i + 1), blkSize * (i + 1), c[i], eb, path[i], &*comSize);
            // std::cout << path[i] << std::endl;
        }
        delete [] ary[i];
    }
    delete [] ary;

	// std::cout << "max:   " << maxBlk<< std::endl;
	int count[maxBlk];
	for(size_t i = 0; i < maxBlk; ++i) {
		count[i] = 0;
	}
	for (int z = blkNum - 1; z > -1; z--){
		for (int y = blkNum - 1; y > -1; y--){
			for (int x = blkNum - 1; x > -1; x--){
				size_t size = take[x + (y * blkNum) + (z * blkNum * blkNum)];
				if (size > 0)
				{
					for (size_t k = (z - size + 1) * blkSize; k < (z + 1) * blkSize; k++){
						for (size_t j = (y - size + 1) * blkSize; j < (y + 1) * blkSize; j++){
							for (size_t i = (x - size + 1) * blkSize; i < (x + 1)  * blkSize; i++){
								densgrid[i + (j * grid) + (k * grid * grid)] = ts[size - 1][count[size - 1]];
								count[size - 1] ++ ;
							}
						}
					}
				}
			}
		}
	}
	return densgrid;
}

double *Gsp(double *densgrid, char* tp3d, size_t grid_x, size_t grid_y, size_t grid_z,size_t blkSize, size_t meat, size_t bread, char* inName, double eb, char* cfg, size_t lvl, size_t* comSize){
    size_t blkNum_x = grid_x / blkSize;
	size_t blkNum_y = grid_y / blkSize;
	size_t blkNum_z = grid_z / blkSize;
    char * merge = new char[grid_x*grid_y*grid_z]{0};

    for (size_t z = 0; z < blkNum_z; z++){
		for (size_t y = 0; y < blkNum_y; y++){
			for (size_t x = 0; x < blkNum_x; x++){
				if (tp3d[x + (y * blkNum_x) + (z * blkNum_x * blkNum_y)] == 0){
					double avg;
					if (x > 0 && tp3d[(x - 1) + (y * blkNum_x) + (z * blkNum_x * blkNum_y)] == 1)
					{
						for (size_t k = z * blkSize; k < (z + 1) * blkSize ; k++){
							for (size_t j = y * blkSize; j < (y + 1) * blkSize ; j++){
								for (size_t p = 1; p <= meat ; p++){
									avg += densgrid[(x * blkSize) - p + (j * grid_x) + (k * grid_x * grid_y)];	
								}
								for (size_t i = x * blkSize; i < x * blkSize + bread; i++){
									size_t cpos = i + (j * grid_x) + (k * grid_x * grid_y);
									if (merge[cpos] == 0){
										densgrid[cpos] = avg / meat;
										// densgrid[cpos] = 2;
										merge[cpos] += 1;
									}
									else if (merge[cpos] == 1){
										densgrid[cpos] = (densgrid[cpos] + avg/meat) / 2;
										// densgrid[cpos] = 4;
										merge[cpos] += 1;
									}else{
										densgrid[cpos] = (densgrid[cpos] * 2 + avg/meat) / 3;
										// densgrid[cpos] = 8;
									}
								}
								avg = 0;			
							}
						}			
					}
					if (x < blkNum_x - 1 && tp3d[(x + 1) + (y * blkNum_x) + (z * blkNum_x * blkNum_y)] == 1)
					{
						for (size_t k = z * blkSize; k < (z + 1) * blkSize ; k++){
							for (size_t j = y * blkSize; j < (y + 1) * blkSize ; j++){
								for (size_t p = 1; p <= meat ; p++){
									avg += densgrid[((x + 1) * blkSize) + p - 1 + (j * grid_x) + (k * grid_x * grid_y)];	
								}
								for (size_t i = (x + 1) * blkSize - bread; i < (x + 1) * blkSize ; i++){
									size_t cpos = i + (j * grid_x) + (k * grid_x * grid_y);
									if (merge[cpos] == 0){
										densgrid[cpos] = avg / meat;
										// densgrid[cpos] = 2;
										merge[cpos] += 1;
									}
									else if (merge[cpos] == 1){
										densgrid[cpos] = (densgrid[cpos] + avg/meat) / 2;
										// densgrid[cpos] = 4;
										merge[cpos] += 1;
									}else{
										densgrid[cpos] = (densgrid[cpos] * 2 + avg/meat) / 3;
										// densgrid[cpos] = 8;
									}
								}
								avg = 0;			
							}
						}			
					}
					if (y > 0 && tp3d[x  + ((y - 1) * blkNum_x) + (z * blkNum_x * blkNum_y)] == 1)
					{
						for (size_t k = z * blkSize; k < (z + 1) * blkSize ; k++){
							for (size_t i = x * blkSize; i < (x + 1) * blkSize ; i++){
								for (size_t p = 1; p <= meat; p++){
									avg += densgrid[i + ((y * blkSize - p) * grid_x) + (k * grid_x * grid_y)];	
								}
								for (size_t j = y * blkSize; j < y * blkSize + bread; j++){
									size_t cpos = i + (j * grid_x) + (k * grid_x * grid_y);
									if (merge[cpos] == 0){
										densgrid[cpos] = avg / meat;
										// densgrid[cpos] = 2;
										merge[cpos] += 1;
									}
									else if (merge[cpos] == 1){
										densgrid[cpos] = (densgrid[cpos] + avg/meat) / 2;
										// densgrid[cpos] = 4;
										merge[cpos] += 1;
									}else{
										densgrid[cpos] = (densgrid[cpos] * 2 + avg/meat) / 3;
										// densgrid[cpos] = 8;
									}
								}
								avg = 0;			
							}
						}			
					}
					if (y < blkNum_y - 1 && tp3d[x  + ((y + 1) * blkNum_x) + (z * blkNum_x * blkNum_y)] == 1)
					{
						for (size_t k = z * blkSize; k < (z + 1) * blkSize ; k++){
							for (size_t i = x * blkSize; i < (x + 1) * blkSize ; i++){
								for (size_t p = 1; p <= meat ; p++){
									avg += densgrid[i + ((y + 1) * blkSize + p - 1) * grid_x + (k * grid_x * grid_y)];	
								}
								for (size_t j = (y + 1) * blkSize - bread; j < (y + 1) * blkSize; j++){
									size_t cpos = i + (j * grid_x) + (k * grid_x * grid_y);
									if (merge[cpos] == 0){
										densgrid[cpos] = avg / meat;
										// densgrid[cpos] = 2;
										merge[cpos] += 1;
									}
									else if (merge[cpos] == 1){
										densgrid[cpos] = (densgrid[cpos] + avg/meat) / 2;
										// densgrid[cpos] = 4;
										merge[cpos] += 1;
									}else{
										densgrid[cpos] = (densgrid[cpos] * 2 + avg/meat) / 3;
										// densgrid[cpos] = 8;
									}
								}
								avg = 0;			
							}
						}			
					}
					if (z > 0 && tp3d[x  + (y * blkNum_x) + ((z - 1) * blkNum_x * blkNum_y)] == 1)
					{
						for (size_t j = y * blkSize; j < (y + 1) * blkSize ; j++){
							for (size_t i = x * blkSize; i < (x + 1) * blkSize ; i++){
								for (size_t p = 1; p <= meat; p++){
									avg += densgrid[i + (j * grid_x) + ((z * blkSize - p) * grid_x * grid_y)];	
								}
								for (size_t k = z * blkSize; k < z * blkSize + bread; k++){
									size_t cpos = i + (j * grid_x) + (k * grid_x * grid_y);
									if (merge[cpos] == 0){
										densgrid[cpos] = avg / meat;
										// densgrid[cpos] = 2;
										merge[cpos] += 1;
									}
									else if (merge[cpos] == 1){
										densgrid[cpos] = (densgrid[cpos] + avg/meat) / 2;
										// densgrid[cpos] = 4;
										merge[cpos] += 1;
									}else{
										densgrid[cpos] = (densgrid[cpos] * 2 + avg/meat) / 3;
										// densgrid[cpos] = 8;
									}
								}
								avg = 0;			
							}
						}			
					}
					if (z < blkNum_z - 1 && tp3d[x  + (y * blkNum_x) + ((z + 1) * blkNum_x * blkNum_y)] == 1)
					{
						for (size_t j = y * blkSize; j < (y + 1) * blkSize ; j++){
							for (size_t i = x * blkSize; i < (x + 1) * blkSize ; i++){
								for (size_t p = 1; p <= meat; p++){
									avg += densgrid[i + (j * grid_x) + ((z + 1) * blkSize + p - 1) * grid_x * grid_y];	
								}
								for (size_t k = (z + 1) * blkSize - bread; k < (z + 1) * blkSize; k++){
									size_t cpos = i + (j * grid_x) + (k * grid_x * grid_y);
									if (merge[cpos] == 0){
										densgrid[cpos] = avg / meat;
										// densgrid[cpos] = 2;
										merge[cpos] += 1;
									}
									else if (merge[cpos] == 1){
										densgrid[cpos] = (densgrid[cpos] + avg/meat) / 2;
										// densgrid[cpos] = 4;
										merge[cpos] += 1;
									}else{
										densgrid[cpos] = (densgrid[cpos] * 2 + avg/meat) / 3;
										// densgrid[cpos] = 8;
									}
								}
								avg = 0;			
							}
						}			
					}
				}
            }
		}	
	}	
	delete[] merge;
	char path[650];
	sprintf(path, "%s.3d_%d", inName, lvl);
    densgrid = AMR_compress(densgrid, cfg, grid_x, grid_y, grid_z, 0, eb, path, &*comSize);
    return densgrid;
}

struct Node
{
    size_t count; 
	size_t xBegin;
	size_t yBegin;
	size_t zBegin;
	size_t xEnd;
	size_t yEnd;
	size_t zEnd;
    Node *left, *right;
};
  
Node *insertRec(char* arr, size_t xBegin, size_t xEnd, size_t yBegin, size_t yEnd, size_t zBegin, size_t zEnd, size_t count, size_t a, size_t b, size_t c, size_t d, size_t length)
{
    struct Node* node = new Node;
    node->count = count;
	node->xBegin = xBegin;
	node->yBegin = yBegin;
	node->zBegin = zBegin;
	node->xEnd = xEnd;
	node->yEnd = yEnd;
	node->zEnd = zEnd;
    
    // std::cout << a << " " << b << " " << c << " " << d << std::endl;
    size_t xSize = xEnd - xBegin + 1;
    size_t ySize = yEnd - yBegin + 1;
    size_t zSize = zEnd - zBegin + 1;
    size_t size = xSize * ySize * zSize;
	
	if (count == size || count == 0){

    // if ((double)count/size >= 0.95 || count == 0){
	// 	if (count != 0)
	// 	{
	// 		node->count = size;
	// 	}
    } else {
        if (xSize == ySize && ySize == zSize){
            size_t* blkCount = new size_t[8]{0};
            for (size_t z = 0; z <= 1; z++){
                for (size_t y = 0; y <= 1; y++){
                    for (size_t x = 0; x <= 1; x++){
                        for (size_t k = zBegin + z * zSize/2; k < zBegin + (z + 1) * zSize/2; k++){
                            for (size_t j = yBegin + y * ySize/2; j < yBegin + (y + 1) * ySize/2; j++){
                                for (size_t i = xBegin + x * xSize/2; i < xBegin + (x + 1) * xSize/2; i++){
									// std::cout << i << " " << j << " " << k << std::endl;
                                    if (arr[i + (j * length) + (k * length * length)] == 1){
                                        blkCount[x + (y * 2) + (z * 2 * 2)]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // for (size_t i = 0; i < 8; i++)
            // {
            //     std::cout << blkCount[i] << " ";
            // }
            size_t xRate = abs(int((blkCount[0] + blkCount[2] + blkCount[4] + blkCount[6]) - (blkCount[1] + blkCount[3] + blkCount[5] + blkCount[7])));
            size_t yRate = abs(int((blkCount[0] + blkCount[1] + blkCount[4] + blkCount[5]) - (blkCount[2] + blkCount[3] + blkCount[6] + blkCount[7])));
            size_t zRate = abs(int((blkCount[0] + blkCount[1] + blkCount[2] + blkCount[3]) - (blkCount[4] + blkCount[5] + blkCount[6] + blkCount[7])));
            size_t maxRate = max(max(xRate, yRate), zRate);
            // std::cout << maxRate << std::endl;

            if (maxRate == xRate)
            {
                //std::cout << "x cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xBegin + xSize/2 - 1, yBegin, yEnd, zBegin, zEnd, blkCount[0] + blkCount[2] + blkCount[4] + blkCount[6], blkCount[0],  blkCount[2], blkCount[4], blkCount[6], length);
                node->right = insertRec(arr, xBegin + xSize/2, xEnd, yBegin, yEnd, zBegin, zEnd, blkCount[1] + blkCount[3] + blkCount[5] + blkCount[7], blkCount[1],  blkCount[3], blkCount[5], blkCount[7], length);
            } else if (maxRate == yRate) {
                //std::cout << "y cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xEnd, yBegin, yBegin + ySize/2 - 1, zBegin, zEnd, blkCount[0] + blkCount[1] + blkCount[4] + blkCount[5], blkCount[0],  blkCount[1], blkCount[4], blkCount[5], length);
                node->right = insertRec(arr, xBegin, xEnd, yBegin + ySize/2, yEnd, zBegin, zEnd, blkCount[2] + blkCount[3] + blkCount[6] + blkCount[7], blkCount[2],  blkCount[3], blkCount[6], blkCount[7], length);
            } else {
                //std::cout << "z cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin, zBegin + zSize/2 - 1, blkCount[0] + blkCount[1] + blkCount[2] + blkCount[3], blkCount[0],  blkCount[1], blkCount[2], blkCount[3], length);
                node->right = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin + zSize/2, zEnd, blkCount[4] + blkCount[5] + blkCount[6] + blkCount[7], blkCount[4],  blkCount[5], blkCount[6], blkCount[7], length);
            }
        } else if (xSize == ySize && ySize > zSize) {
            size_t* blkCount = new size_t[4]{a, b, c, d};

            size_t xRate = abs(int((blkCount[0] + blkCount[2]) - (blkCount[1] + blkCount[3])));
            size_t yRate = abs(int((blkCount[0] + blkCount[1]) - (blkCount[2] + blkCount[3])));
            size_t maxRate = max(xRate, yRate);
            if (maxRate == xRate)
            {
                //std::cout << "x cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xBegin + xSize/2 - 1, yBegin, yEnd, zBegin, zEnd, blkCount[0] + blkCount[2], blkCount[0], blkCount[2], 0, 0, length);
                node->right = insertRec(arr, xBegin + xSize/2, xEnd, yBegin, yEnd, zBegin, zEnd, blkCount[1] + blkCount[3], blkCount[1], blkCount[3], -1, -1, length);
            } else {
                //std::cout << "y cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xEnd, yBegin, yBegin + ySize/2 - 1, zBegin, zEnd, blkCount[0] + blkCount[1], blkCount[0], blkCount[1], 0, 0, length);
                node->right = insertRec(arr, xBegin, xEnd, yBegin + ySize/2, yEnd, zBegin, zEnd, blkCount[2] + blkCount[3], blkCount[2], blkCount[3], 0, 0, length);
            }
        }
        else if (xSize == zSize && xSize > ySize) {
            size_t* blkCount = new size_t[4]{a, b, c, d};

            size_t xRate = abs(int((blkCount[0] + blkCount[2]) - (blkCount[1] + blkCount[3])));
            size_t zRate = abs(int((blkCount[0] + blkCount[1]) - (blkCount[2] + blkCount[3])));
            size_t maxRate = max(xRate, zRate);
            if (maxRate == xRate)
            {
                //std::cout << "x cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xBegin + xSize/2 - 1, yBegin, yEnd, zBegin, zEnd, blkCount[0] + blkCount[2], blkCount[0], blkCount[2], 0, 0, length);
                node->right = insertRec(arr, xBegin + xSize/2, xEnd, yBegin, yEnd, zBegin, zEnd, blkCount[1] + blkCount[3], blkCount[1], blkCount[3], -1, -1, length);
            } else {
                //std::cout << "z cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin, zBegin + zSize/2 - 1, blkCount[0] + blkCount[1], blkCount[0], blkCount[1], 0, 0, length);
                node->right = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin + zSize/2, zEnd, blkCount[2] + blkCount[3], blkCount[2], blkCount[3], 0, 0, length);
            }
        }
        else if (ySize == zSize && ySize > xSize) {
            size_t* blkCount = new size_t[4]{a, b, c, d};

            size_t yRate = abs(int((blkCount[0] + blkCount[2]) - (blkCount[1] + blkCount[3])));
            size_t zRate = abs(int((blkCount[0] + blkCount[1]) - (blkCount[2] + blkCount[3])));
            size_t maxRate = max(yRate, zRate);
            if (maxRate == yRate)
            {
                //std::cout << "y cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xEnd, yBegin, yBegin + ySize/2 - 1, zBegin, zEnd, blkCount[0] + blkCount[2], blkCount[0], blkCount[2], 0, 0, length);
                node->right = insertRec(arr, xBegin, xEnd, yBegin + ySize/2, yEnd, zBegin, zEnd, blkCount[1] + blkCount[3], blkCount[1], blkCount[3], -1, -1, length);
            } else {
                //std::cout << "z cut"<< std::endl;
                node->left  = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin, zBegin + zSize/2 - 1, blkCount[0] + blkCount[1], blkCount[0], blkCount[1], 0, 0, length);
                node->right = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin + zSize/2, zEnd, blkCount[2] + blkCount[3], blkCount[2], blkCount[3], 0, 0, length);
            }
        }
        else if (xSize == ySize && ySize < zSize) {
            size_t* blkCount = new size_t[2]{a, b};
            node->left  = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin, zBegin + zSize/2 - 1, blkCount[0], 0, 0, 0, 0, length);
            node->right = insertRec(arr, xBegin, xEnd, yBegin, yEnd, zBegin + zSize/2, zEnd, blkCount[1], 0, 0, 0, 0, length);
        }
        else if (xSize == zSize && zSize < ySize) {
            size_t* blkCount = new size_t[2]{a, b};
            node->left  = insertRec(arr, xBegin, xEnd, yBegin, yBegin + ySize/2 - 1, zBegin, zEnd, blkCount[0], 0, 0, 0, 0, length);
            node->right = insertRec(arr, xBegin, xEnd, yBegin + ySize/2, yEnd, zBegin, zEnd, blkCount[1], 0, 0, 0, 0, length);
        }
        else if (ySize == zSize && ySize < xSize) {
            size_t* blkCount = new size_t[2]{a, b};
            node->left  = insertRec(arr, xBegin, xBegin + xSize/2 - 1, yBegin, yEnd, zBegin, zEnd, blkCount[0], 0, 0, 0, 0, length);
            node->right = insertRec(arr, xBegin + xSize/2, xEnd, yBegin, yEnd, zBegin, zEnd, blkCount[1], 0, 0, 0, 0, length);
        }
    }
    return node;
}

void countLeafNodes(Node *root, size_t* arr)
{
    // if node is null, return
    if (!root)
        return;

    if (!root->left && !root->right && root->count > 0)
    {
		size_t dep = log2(root->count);
        arr[dep] += 1;
        return;
    }

    // if child exists, check for leaf recursively
    if (root->left)
       countLeafNodes(root->left, arr);
         
    if (root->right)
       countLeafNodes(root->right, arr);
}

void tree_2_4d(Node *root, double* data, double** arr, size_t* leafCnt, size_t* ind, size_t grid, size_t blkSize, char*** reTree, size_t* bigInd)
{
    // if node is null, return
    if (!root)
        return;

    if (!root->left && !root->right && root->count > 0)
    {
		
		size_t dep = log2(root->count);

		reTree[dep][bigInd[dep]][0] = root->xBegin;
		reTree[dep][bigInd[dep]][1] = root->xEnd;
		reTree[dep][bigInd[dep]][2] = root->yBegin;
		reTree[dep][bigInd[dep]][3] = root->yEnd;
		reTree[dep][bigInd[dep]][4] = root->zBegin;
		reTree[dep][bigInd[dep]][5] = root->zEnd;
		bigInd[dep]++;
		// std::cout << "qqqqqqqqqqq" << std::endl;

	    size_t xSize = root->xEnd - root->xBegin + 1;
        size_t ySize = root->yEnd - root->yBegin + 1;
        size_t zSize = root->zEnd - root->zBegin + 1;
		// for (size_t z = root->zBegin * blkSize; z < (root->zEnd + 1) * blkSize; z++){
		// 	for (size_t y = root->yBegin * blkSize; y < (root->yEnd + 1) * blkSize; y++){
		// 		for (size_t x = root->xBegin * blkSize; x < (root->xEnd + 1) * blkSize; x++){
		// 			// std::cout<< root->count << " " <<dep << std::endl;
		// 			arr[dep][ind[dep]] = data[x + (y * grid) + (z * grid * grid)];
		// 			ind[dep]++;
		// 			// std::cout<< arr[dep][ind[dep]] << std::endl;
		// 		}
		// 	}
		// }
		// std::cout<< dep << std::endl;
		
		if (xSize == ySize)
		{
			for (size_t z = root->zBegin * blkSize; z < (root->zEnd + 1) * blkSize; z++){
				for (size_t y = root->yBegin * blkSize; y < (root->yEnd + 1) * blkSize; y++){
					for (size_t x = root->xBegin * blkSize; x < (root->xEnd + 1) * blkSize; x++){
						// std::cout<< root->count << " " <<dep << std::endl;
						arr[dep][ind[dep]] = data[x + (y * grid) + (z * grid * grid)];
						ind[dep]++;
						// std::cout<< arr[dep][ind[dep]] << std::endl;
					}
				}
			}
		} else if (xSize == zSize){
			for (size_t y = root->yBegin * blkSize; y < (root->yEnd + 1) * blkSize; y++){
				for (size_t z = root->zBegin * blkSize; z < (root->zEnd + 1) * blkSize; z++){
					for (size_t x = root->xBegin * blkSize; x < (root->xEnd + 1) * blkSize; x++){
						// std::cout<< root->count << " " <<dep << std::endl;
						arr[dep][ind[dep]] = data[x + (y * grid) + (z * grid * grid)];
						ind[dep]++;
						// std::cout<< arr[dep][ind[dep]] << std::endl;
					}
				}
			}
		} else {
			for (size_t x = root->xBegin * blkSize; x < (root->xEnd + 1) * blkSize; x++){
				for (size_t z = root->zBegin * blkSize; z < (root->zEnd + 1) * blkSize; z++){
					for (size_t y = root->yBegin * blkSize; y < (root->yEnd + 1) * blkSize; y++){
						// std::cout<< root->count << " " <<dep << std::endl;
						arr[dep][ind[dep]] = data[x + (y * grid) + (z * grid * grid)];
						ind[dep]++;
						// std::cout<< arr[dep][ind[dep]] << std::endl;
					}
				}
			}
		}
        return;
    }

    // if child exists, check for leaf recursively
    if (root->left)
       tree_2_4d(root->left, data, arr, leafCnt, ind, grid, blkSize, reTree, bigInd);
         
    if (root->right)
       tree_2_4d(root->right, data, arr, leafCnt, ind, grid, blkSize, reTree, bigInd);
}

void treeBack(double* densgrid, double** tree, size_t grid, size_t blkSize, char*** reTree, size_t* leafCnt, char maxTreeBlk)
{
	int count[maxTreeBlk];
	for(size_t i = 0; i < maxTreeBlk; ++i) {
		count[i] = 0;
	}

	size_t xBegin, xEnd, yBegin, yEnd, zBegin, zEnd, xSize, ySize, zSize; 
	for (size_t i = 0; i < maxTreeBlk; i++){
		for (size_t j = 0; j < leafCnt[i]; j++){
			xBegin = reTree[i][j][0];
			xEnd   = reTree[i][j][1];
			yBegin = reTree[i][j][2];
			yEnd   = reTree[i][j][3];
			zBegin = reTree[i][j][4];
			zEnd   = reTree[i][j][5];
			xSize = xEnd - xBegin + 1;
            ySize = yEnd - yBegin + 1;
            zSize = zEnd - zBegin + 1;
			// for (size_t z = zBegin * blkSize; z < (zEnd+ 1) * blkSize; z++){
			// 	for (size_t y = yBegin * blkSize; y < (yEnd + 1) * blkSize; y++){
			// 		for (size_t x = xBegin * blkSize; x < (xEnd + 1) * blkSize; x++){
			// 			densgrid[(x + (y * grid) + (z * grid * grid))] = tree[i][count[i]];
			// 			count[i]++;
			// 		}
			// 	}
			// }
			if (xSize == ySize)
			{
				for (size_t z = zBegin * blkSize; z < (zEnd+ 1) * blkSize; z++){
					for (size_t y = yBegin * blkSize; y < (yEnd + 1) * blkSize; y++){
						for (size_t x = xBegin * blkSize; x < (xEnd + 1) * blkSize; x++){
							densgrid[(x + (y * grid) + (z * grid * grid))] = tree[i][count[i]];
							count[i]++;
						}
					}
				}
			} else if (xSize == zSize) {
				for (size_t y = yBegin * blkSize; y < (yEnd + 1) * blkSize; y++){
					for (size_t z = zBegin * blkSize; z < (zEnd+ 1) * blkSize; z++){
						for (size_t x = xBegin * blkSize; x < (xEnd + 1) * blkSize; x++){
							densgrid[(x + (y * grid) + (z * grid * grid))] = tree[i][count[i]];
							count[i]++;
						}
					}
				}
			} else {
				for (size_t x = xBegin * blkSize; x < (xEnd + 1) * blkSize; x++){
					for (size_t z = zBegin * blkSize; z < (zEnd+ 1) * blkSize; z++){
						for (size_t y = yBegin * blkSize; y < (yEnd + 1) * blkSize; y++){
							densgrid[(x + (y * grid) + (z * grid * grid))] = tree[i][count[i]];
							count[i]++;
						}
					}
				}
			}
		}
	}
	
	return;
}

double *ADKtree (double *densgrid, size_t cnt, char* tp3d, size_t grid, size_t blkNum, size_t blkSize, char* inName, double eb, char* cfg, size_t lvl, size_t* comSize) {
	size_t count = cnt;
	struct Node *root = NULL;
	root = insertRec(tp3d, 0, blkNum - 1, 0, blkNum - 1, 0, blkNum - 1, count, 0, 0, 0, 0, blkNum);

	size_t maxDep = log2(blkNum*blkNum*blkNum)+1;
	size_t* leafCnt = new size_t[maxDep]{0};
	countLeafNodes(root, leafCnt);

	size_t maxTreeBlk = maxDep;
	double **tree = new double*[maxTreeBlk];
    for(int i = 0; i < maxTreeBlk; ++i) {
		size_t treeBlkSize = pow(2, i);
        tree[i] = new double[leafCnt[i] * treeBlkSize * blkSize * blkSize * blkSize];
    }

	size_t *ind = new size_t[maxTreeBlk]{0};
	size_t *bigInd = new size_t[maxTreeBlk]{0};

	char*** reTree = new char**[maxTreeBlk];
	for(size_t i = 0; i < maxTreeBlk; ++i) {
        reTree[i] = new char* [leafCnt[i]];
		for (size_t j = 0; j < leafCnt[i]; j++)
		{
			reTree[i][j] = new char [6];
		}		
    }

	tree_2_4d(root, densgrid, tree, leafCnt, ind, grid, blkSize, reTree, bigInd);

	char path[maxTreeBlk][650];
	for(size_t i = 0; i < maxTreeBlk; ++i) {
		if (leafCnt[i] > 0){
			sprintf(path[i], "%s.3d_%d_%d", inName, lvl , i);
			double treeBlkSize = pow(2, i);
			double temp_0 = std::cbrt(treeBlkSize);
			double temp_1 = std::cbrt(treeBlkSize * 2);
			double temp_2 = std::cbrt(treeBlkSize/2);
			size_t h = 0;
			if (temp_0 - int(temp_0) == 0){
				size_t blkNum_x = int(temp_0);
				if (leafCnt[i] == 1){
					h = 0;
				} else {
					h = leafCnt[i];
				}
				tree[i] = AMR_compress(tree[i], cfg, blkSize * blkNum_x, blkSize * blkNum_x, blkSize * blkNum_x, h, eb, path[i], &*comSize);
			}

			if (temp_1 - int(temp_1) == 0){
				size_t blkNum_x = int(temp_1);
				size_t blkNum_z = blkNum_x/2;
				if (leafCnt[i] == 1){
					h = 0;
				} else {
					h = leafCnt[i];
				}
				tree[i] = AMR_compress(tree[i], cfg, blkSize * blkNum_x, blkSize * blkNum_x, blkSize * blkNum_x/2, h, eb, path[i], &*comSize);
			}

			if (temp_2 - int(temp_2) == 0){
				size_t blkNum_x = int(temp_2);
				size_t blkNum_z = blkNum_x * 2;

				if (leafCnt[i] == 1){
					h = 0;
				} else {
					h = leafCnt[i];
				}
				tree[i] = AMR_compress(tree[i], cfg, blkSize * blkNum_x, blkSize * blkNum_x, blkSize * 2 * blkNum_x, h, eb, path[i], &*comSize);
			}
		}
    }

	treeBack(densgrid, tree, grid, blkSize, reTree, leafCnt,maxTreeBlk);
	return densgrid;
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
	// num_cells = 900000;
	std::cout << "Number of cells in the file: " << num_cells << std::endl;

	double* xx = new double[num_cells];
	double* yy = new double[num_cells];
	double* zz = new double[num_cells];
    double* oriVol = new double[num_cells];
	double* vol = new double[num_cells];
	double* rho = new double[num_cells];
	double* temp = new double[num_cells];

    size_t size[4] = {0, 0, 0, 0};

    double minRho =FLT_MAX;
	double max[5] = { -FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX };

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

		if (rho[i] > max[0])
			max[0] = rho[i];
		if (rho[i] < minRho)
			minRho = rho[i];
		
		if (xx[i] > max[1])
			max[1] = xx[i];

		if (yy[i] > max[2])
			max[2] = yy[i];

		if (zz[i] > max[3])
			max[3] = zz[i];

		vol[i] = std::cbrt(oriVol[i]);

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
	
	double lVol[iden];

	std::cout << "number of AMR levels: " << iden << std::endl;

	for (size_t i = 0; i < iden; i++)
	{
		lVol[i] = lVoltemp[i];
	}
	std::sort(lVol, lVol + iden);

	size_t rate[4] = {1, 0, 0, 0};

	for (size_t i = 1; i < iden; i++)
	{
		rate[i] = lVol[i]/lVol[0];
	}

	size_t grid[iden];
	for (size_t i = 0; i < iden; i++)
	{
		grid[i] = std::max(std::max(max[1], max[2]), max[3]);
		grid[i] = grid[i]/lVol[i];
		grid[i] = pow(2, ceil(log2(grid[i])));
	}

	std::cout << "output dims: " << " ";
	for (int i = 0; i < iden; ++i)
        std::cout << grid[i] << " ";
	std::cout << std::endl;
	
	size_t blkSize[4] = {16, 8, 4, 2};
	size_t blkNum = grid[0] / blkSize[0];

	size_t pad[iden];
	pad[0] = 4;
	for (int i = 1; i < iden; ++i)
        pad[i] = blkSize[i]/2;


	double range = max[0] - minRho;
	std::cout << "range:" << range << std::endl;

	double **densgrid = new double*[iden];
    for(int i = 0; i < iden; ++i) {
        densgrid[i] = new double[grid[i]*grid[i]*grid[i]]{0.0};
    }

	double **origrid = new double*[iden];
    for(int i = 0; i < iden; ++i) {
        origrid[i] = new double[grid[i]*grid[i]*grid[i]]{0.0};
    }

	char **tp3d = new char*[iden];
    for(int i = 0; i < iden; ++i) {
        tp3d[i] = new char[blkNum*blkNum*blkNum]{0};
    }

	//to 3d
	for (size_t i = 0; i < num_cells; ++i)
	{
		size_t xi, yi, zi, cpos;
		// Debug: Make sure 'round' is working right
		double dxi = (xx[i] / lVol[0]);
		if (round(dxi) != dxi)
		{
			std::cout.precision(10);
			if( abs(dxi - round(dxi)) > 0.5)
				cout << "WARNING: Potential mapping error: " << dxi << " vs " << round(dxi) << std::endl;
			
		}

		xi = round(xx[i] / lVol[0]);
		yi = round(yy[i] / lVol[0]);
		zi = round(zz[i] / lVol[0]);

		if (vol[i] == lVol[0])
		{
			cpos = xi + (yi * (grid[0])) + (zi * (grid[0]) * (grid[0]));
			densgrid[0][cpos] = rho[i];
			size[0]++;
			origrid[0][cpos] = rho[i];
		}
		else if (vol[i] == lVol[1])
		{
            cpos = xi/rate[1] + (yi/rate[1] * (grid[1])) + (zi/rate[1] * (grid[1]) * (grid[1]));
            densgrid[1][cpos] = rho[i];
			size[1]++;
			origrid[1][cpos] = rho[i];
		}
		else if (iden > 2 && vol[i] == lVol[2])
		{
            cpos = xi/rate[2] + (yi/rate[2] * (grid[2])) + (zi/rate[2] * (grid[2]) * (grid[2]));
            densgrid[2][cpos] = rho[i];
			size[2]++;
			origrid[2][cpos] = rho[i];
		}
		else
		{
            cpos = xi/rate[3] + (yi/rate[3] * (grid[3])) + (zi/rate[3] * (grid[3]) * (grid[3]));
            densgrid[3][cpos] = rho[i];
			size[3]++;
			origrid[3][cpos] = rho[i];
		}
	}
	
	size_t cnt[iden];
	for (size_t i = 0; i < iden; i++)
	{
		cnt[i] = 0;
	}
	//count
    if (iden == 2)
	{
		for (size_t z = 0; z < blkNum; z++){
			for (size_t y = 0; y < blkNum; y++){
				for (size_t x = 0; x < blkNum; x++){
					size_t ind_1 = 0;
					for (size_t k = z * blkSize[1]; k < (z + 1) * blkSize[1]; k++){
						for (size_t j = y * blkSize[1]; j < (y + 1) * blkSize[1]; j++){
							for (size_t i = x * blkSize[1]; i < (x + 1) * blkSize[1]; i++){
								if (densgrid[1][i + (j * grid[1]) + (k * grid[1] * grid[1])] != 0){
									ind_1 = 1;
								}
							}
						}
					}
					if (ind_1 == 0){
						tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] = 0;
						tp3d[0][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[0]++;
					}else{
						tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[1]++;
					}
				}
			}	
		}
		std::cout << "lvl 0 unit block number: " << cnt[0]<< std::endl;
		std::cout << "lvl 1 unit block number: " << cnt[1]<< std::endl;
	} else if (iden == 3) {
		for (size_t z = 0; z < blkNum; z++){
			for (size_t y = 0; y < blkNum; y++){
				for (size_t x = 0; x < blkNum; x++){
					size_t ind_1 = 0;
					for (size_t k = z * blkSize[1]; k < (z + 1) * blkSize[1]; k++){
						for (size_t j = y * blkSize[1]; j < (y + 1) * blkSize[1]; j++){
							for (size_t i = x * blkSize[1]; i < (x + 1) * blkSize[1]; i++){
								if (densgrid[1][i + (j * grid[1]) + (k * grid[1] * grid[1])] != 0){
									ind_1 = 1;
								}
							}
						}
					}
					if (ind_1 == 0){
						tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] = 0;
					}else{
						tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[1]++;
					}

					size_t ind_2 = 0;
					for (size_t k = z * blkSize[2]; k < (z + 1) * blkSize[2]; k++){
						for (size_t j = y * blkSize[2]; j < (y + 1) * blkSize[2]; j++){
							for (size_t i = x * blkSize[2]; i < (x + 1) * blkSize[2]; i++){
								if (densgrid[2][i + (j * grid[2]) + (k * grid[2] * grid[2])] != 0){
									ind_2 = 1;
								}
							}
						}
					}
					if (ind_2 == 0){
						tp3d[2][x + (y * blkNum) + (z * blkNum * blkNum)] = 0;
					}else{
						tp3d[2][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[2]++;
					}

					if (tp3d[2][x + (y * blkNum) + (z * blkNum * blkNum)] == 0 && tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] == 0){
						tp3d[0][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[0]++;
					}
				}
			}	
		}
		std::cout << "lvl 0 unit block number: " << cnt[0]<< std::endl;
		std::cout << "lvl 1 unit block number: " << cnt[1]<< std::endl;
		std::cout << "lvl 2 unit block number: " << cnt[2]<< std::endl;
	}  else  {
		for (size_t z = 0; z < blkNum; z++){
			for (size_t y = 0; y < blkNum; y++){
				for (size_t x = 0; x < blkNum; x++){
					size_t ind_1 = 0;
					for (size_t k = z * blkSize[1]; k < (z + 1) * blkSize[1]; k++){
						for (size_t j = y * blkSize[1]; j < (y + 1) * blkSize[1]; j++){
							for (size_t i = x * blkSize[1]; i < (x + 1) * blkSize[1]; i++){
								if (densgrid[1][i + (j * grid[1]) + (k * grid[1] * grid[1])] != 0){
									ind_1 = 1;
								}
							}
						}
					}
					if (ind_1 == 0){
						tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] = 0;
					}else{
						tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[1]++;
					}

					size_t ind_2 = 0;
					for (size_t k = z * blkSize[2]; k < (z + 1) * blkSize[2]; k++){
						for (size_t j = y * blkSize[2]; j < (y + 1) * blkSize[2]; j++){
							for (size_t i = x * blkSize[2]; i < (x + 1) * blkSize[2]; i++){
								if (densgrid[2][i + (j * grid[2]) + (k * grid[2] * grid[2])] != 0){
									ind_2 = 1;
								}
							}
						}
					}
					if (ind_2 == 0){
						tp3d[2][x + (y * blkNum) + (z * blkNum * blkNum)] = 0;
					}else{
						tp3d[2][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[2]++;
					}

					size_t ind_3 = 0;
					for (size_t k = z * blkSize[3]; k < (z + 1) * blkSize[3]; k++){
						for (size_t j = y * blkSize[3]; j < (y + 1) * blkSize[3]; j++){
							for (size_t i = x * blkSize[3]; i < (x + 1) * blkSize[3]; i++){
								if (densgrid[3][i + (j * grid[3]) + (k * grid[3] * grid[3])] != 0){
									ind_3 = 1;
								}
							}
						}
					}
					if (ind_3 == 0){
						tp3d[3][x + (y * blkNum) + (z * blkNum * blkNum)] = 0;
					}else{
						tp3d[3][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[3]++;
					}

					if (tp3d[2][x + (y * blkNum) + (z * blkNum * blkNum)] == 0 && tp3d[1][x + (y * blkNum) + (z * blkNum * blkNum)] == 0 && tp3d[3][x + (y * blkNum) + (z * blkNum * blkNum)] == 0){
						tp3d[0][x + (y * blkNum) + (z * blkNum * blkNum)] = 1;
						cnt[0]++;
					}
				}
			}	
		}
		std::cout << "lvl 0 unit block number: " << cnt[0]<< std::endl;
		std::cout << "lvl 1 unit block number: " << cnt[1]<< std::endl;
		std::cout << "lvl 2 unit block number: " << cnt[2]<< std::endl;
		std::cout << "lvl 3 unit block number: " << cnt[3]<< std::endl;
	} 

 	char *cfg = argv[2];
	double eb[4] = {atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6])};
	size_t outSize = 0;
	for(int i = 0; i < iden; ++i) {
        if (size[i]/pow(grid[i],3) < 0.5) {
			std::cout << "level " << i << " using: " << "OpST" << std::endl;
			densgrid[i] = Opst_compress(densgrid[i], tp3d[i], grid[i], blkSize[i], inName, eb[i], cfg, i, &outSize);
		}
		else if (size[i]/pow(grid[i],3) > 0.6) {
			std::cout << "level " << i << " using: " << "GSP" << std::endl;
			densgrid[i] = Gsp(densgrid[i], tp3d[i], grid[i], grid[i], grid[i], blkSize[i], pad[i], pad[i], inName, eb[i], cfg, i, &outSize);
		} else {
			std::cout << "level " << i << " using: " << "ADKtree" << std::endl;
			densgrid[i] = ADKtree(densgrid[i], cnt[i], tp3d[i], grid[i], blkNum, blkSize[i], inName, eb[i], cfg, i, &outSize);
			// (double *densgrid, size_t cnt, char* tp3d, size_t grid, size_t blkNum, size_t blkSize, char* inName, double eb, char* cfg, size_t lvl)
		}
    }
	std::cout << "Our compression ratio : " << double(num_cells) * 8 / outSize << std::endl;
	std::cout << "Our bitrate : " << double(outSize) * 8 / num_cells << std::endl;

	double normErr[iden] = {0, 0, 0, 0};
	for (size_t j = 0; j < iden; ++j) {
		for (size_t i = 0; i < pow(grid[j],3); ++i) 
		{
			if (origrid[j][i] != 0)
			{
				normErr[j] += (densgrid[j][i] - origrid[j][i]) * (densgrid[j][i] - origrid[j][i]);
			}
		}
	}
	double mse = 0;
	for (size_t i = 0; i < iden; ++i) {
		mse += pow(rate[i],3) * normErr[i];
	}

	mse = mse/(pow(grid[0],3));

	std::cout << "Our PSNR: " << 20*log10(range)-10*log10(mse) << std::endl;

    return 0;
}
