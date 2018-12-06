# 遥感图像的分块处理

## 一、基本原理

在遥感图像处理中一般采取分块读写的方法，因为遥感图像尺寸巨大（长宽可达30000像素以上），不可能开辟一个大内存把整个图像读进来。分块读取的道理一般大家都懂，但处理上却是有学问的。

在大图像处理中磁盘I/O一般是效率的主要瓶颈。因此如何分块的着眼点应该是如何减少磁盘I/O。一般的图像处理系统采取将块分成256 * 256或者512 * 512的块。实际上这是不利于减少磁盘I/O次数。如下图：

![6_1](/img6/6_1.jpg)

假设上图是按照256 * 256进行分块的，一般而言图像文件在磁盘上是按行存贮的，就是从第一行到最后一行的，那么我们很快就可以看到其中弊端，就是当将数据写入到某一块时，其写入顺序是怎样的呢？首先当然是从块的起始地址写，将块的第一行的数据写入，这时你要问块的第二行数据和第一行数据在连续的？很显然，它们不是连续的。那么你读写时必须先移动文件指针，那么读取一块数据你就需要移动256次文件指针。整幅图像就需要至少移动块数 * 256次文件指针，这样的磁盘I/O次数有点惊人。

因此，一般遥感图像处理的分块策略为如下：

![6_2](/img6/6_2.jpg)

这种分块方法是采取一次读取若干行的方法。这种方法有两个好处：首先降低了程序的逻辑复杂度，每块只是行号变化，块的起始列位置都是0，不需要变化；另一方面，这种分块能大大减少磁盘的I/O次数，这是因为每一块的每行数据在磁盘上都是连续的，因此在读写时只需将文件指针定位到块的起始位置，就能实现整块的读写。毫无疑问，这种分块方法比前一种分块方法磁盘I/O次数要少很多。

## 二、任务描述

使用上周学习的IHS方法，对宽幅遥感图像分别使用两种分块方式进行融合：

第一种方式：使用256 * 256 大小的块；

第二种方式：每块的宽度为图像宽度，高度为256像素。

比较两种分块方式的处理效率。

## 三、实验内容

**第一种方式：使用256 * 256 大小的块：**

```
#include "pch.h"
#include <iostream>
#include <cmath>
#include "./gdal/gdal_priv.h"
#pragma comment(lib, "gdal_i.lib")
using namespace std;

int main()
{
	char* mulPath = (char*)"Mul_large.tif";
	char* panPath = (char*)"Pan_large.tif";
	char* fusPath = (char*)"Fus_large.tif";

	GDALAllRegister();

	// basic parameters
	GDALDataset *poMulDS, *poPanDS, *poFusDS;
	int imgXlen, imgYlen;
	int i = 0, m = 0, n = 0;
	float *bandR, *bandG, *bandB;
	float *bandH, *bandS;
	float *bandP;


	// open datasets
	poMulDS = (GDALDataset*)GDALOpenShared(mulPath, GA_ReadOnly);
	poPanDS = (GDALDataset*)GDALOpenShared(panPath, GA_ReadOnly);
	imgXlen = poMulDS->GetRasterXSize();
	imgYlen = poMulDS->GetRasterYSize();
	poFusDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(
		fusPath, imgXlen, imgYlen, 3, GDT_Byte, NULL);


	int XNum, YNum;//x,y方向上块数
	int sizex;//块宽度
	int sizey;//块高度

	//计算y方向上块数,块宽度,块高度
	sizex = 256;
	sizey = 256;
	YNum = (imgYlen - 1) / sizey + 1;
	XNum = (imgXlen - 1) / sizex + 1;

	// allocating memory
	bandR = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandG = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandB = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandP = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandH = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandS = (float*)CPLMalloc(sizex * sizey * sizeof(float));

	for (m = 0; m < YNum; m++)
	{
		for (n = 0; n < XNum; n++)
		{
			sizex = 256;
			sizey = 256;
			if (n == XNum - 1)
			{
				sizex = (imgXlen - 1) % 256 + 1;
			}

			if (m == YNum - 1)
			{
				sizey = (imgYlen - 1) % 256 + 1;
			}
			//读取当前块数据
			poMulDS->GetRasterBand(1)->RasterIO(GF_Read, n * 256, m * 256, sizex,
				sizey, bandR, sizex, sizey, GDT_Float32, 0, 0);
			poMulDS->GetRasterBand(2)->RasterIO(GF_Read, n * 256, m * 256, sizex,
				sizey, bandG, sizex, sizey, GDT_Float32, 0, 0);
			poMulDS->GetRasterBand(3)->RasterIO(GF_Read, n * 256, m * 256, sizex,
				sizey, bandB, sizex, sizey, GDT_Float32, 0, 0);
			poPanDS->GetRasterBand(1)->RasterIO(GF_Read, n * 256, m * 256, sizex,
				sizey, bandP, sizex, sizey, GDT_Float32, 0, 0);

			for (i = 0; i < sizex * sizey; i++)
			{
				bandH[i] = -sqrt(2.0f) / 6.0f*bandR[i] - sqrt(2.0f) / 6.0f*bandG[i] + sqrt(2.0f) / 3.0f*bandB[i];
				bandS[i] = 1.0f / sqrt(2.0f)*bandR[i] - 1 / sqrt(2.0f)*bandG[i];

				bandR[i] = bandP[i] - 1.0f / sqrt(2.0f)*bandH[i] + 1.0f / sqrt(2.0f)*bandS[i];
				bandG[i] = bandP[i] - 1.0f / sqrt(2.0f)*bandH[i] - 1.0f / sqrt(2.0f)*bandS[i];
				bandB[i] = bandP[i] + sqrt(2.0f)*bandH[i];
			}
			poFusDS->GetRasterBand(1)->RasterIO(GF_Write, n * 256, m * 256, sizex, sizey,
				bandR, sizex, sizey, GDT_Float32, 0, 0);
			poFusDS->GetRasterBand(2)->RasterIO(GF_Write, n * 256, m * 256, sizex, sizey,
				bandG, sizex, sizey, GDT_Float32, 0, 0);
			poFusDS->GetRasterBand(3)->RasterIO(GF_Write, n * 256, m * 256, sizex, sizey,
				bandB, sizex, sizey, GDT_Float32, 0, 0);
		}
	}

	CPLFree(bandR);
	CPLFree(bandG);
	CPLFree(bandB);
	CPLFree(bandH);
	CPLFree(bandS);
	CPLFree(bandP);

	GDALClose(poMulDS);
	GDALClose(poPanDS);
	GDALClose(poFusDS);

	return 0;
}
```

**第二种方式：每块的宽度为图像宽度，高度为256像素：**

```
#include "pch.h"
#include <iostream>
#include <cmath>
#include "./gdal/gdal_priv.h"
#pragma comment(lib, "gdal_i.lib")
using namespace std;

int main()
{
	char* mulPath = (char*)"Mul_large.tif";
	char* panPath = (char*)"Pan_large.tif";
	char* fusPath = (char*)"Fus_large.tif";

	GDALAllRegister();

	// basic parameters
	GDALDataset *poMulDS, *poPanDS, *poFusDS;
	int imgXlen, imgYlen;
	int i = 0, m = 0;
	float *bandR, *bandG, *bandB;
	float *bandH, *bandS;
	float *bandP;


	// open datasets
	poMulDS = (GDALDataset*)GDALOpenShared(mulPath, GA_ReadOnly);
	poPanDS = (GDALDataset*)GDALOpenShared(panPath, GA_ReadOnly);
	imgXlen = poMulDS->GetRasterXSize();
	imgYlen = poMulDS->GetRasterYSize();
	poFusDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(
		fusPath, imgXlen, imgYlen, 3, GDT_Byte, NULL);

	int YNum;//y方向上块数
	int sizex;//块宽度
	int sizey;//块高度

	//计算y方向上块数,块宽度,块高度
	sizex = imgXlen;
	sizey = 256;
	YNum = (imgYlen - 1) / sizey + 1;
	

	// allocating memory
	bandR = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandG = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandB = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandP = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandH = (float*)CPLMalloc(sizex * sizey * sizeof(float));
	bandS = (float*)CPLMalloc(sizex * sizey * sizeof(float));


	for (m = 0; m < YNum; m++)
	{
		if (m == YNum - 1)
		{
			sizey = (imgYlen - 1) % 256 + 1;
		}
		//读取当前块数据
		poMulDS->GetRasterBand(1)->RasterIO(GF_Read, 0, m * 256, sizex,
			sizey, bandR, sizex, sizey, GDT_Float32, 0, 0);
		poMulDS->GetRasterBand(2)->RasterIO(GF_Read, 0, m * 256, sizex,
			sizey, bandG, sizex, sizey, GDT_Float32, 0, 0);
		poMulDS->GetRasterBand(3)->RasterIO(GF_Read, 0, m * 256, sizex,
			sizey, bandB, sizex, sizey, GDT_Float32, 0, 0);
		poPanDS->GetRasterBand(1)->RasterIO(GF_Read, 0, m * 256, sizex,
			sizey, bandP, sizex, sizey, GDT_Float32, 0, 0);

		for (i = 0; i < sizex * sizey; i++)
		{
			bandH[i] = -sqrt(2.0f) / 6.0f*bandR[i] - sqrt(2.0f) / 6.0f*bandG[i] + sqrt(2.0f) / 3.0f*bandB[i];
			bandS[i] = 1.0f / sqrt(2.0f)*bandR[i] - 1 / sqrt(2.0f)*bandG[i];

			bandR[i] = bandP[i] - 1.0f / sqrt(2.0f)*bandH[i] + 1.0f / sqrt(2.0f)*bandS[i];
			bandG[i] = bandP[i] - 1.0f / sqrt(2.0f)*bandH[i] - 1.0f / sqrt(2.0f)*bandS[i];
			bandB[i] = bandP[i] + sqrt(2.0f)*bandH[i];
		}
		poFusDS->GetRasterBand(1)->RasterIO(GF_Write, 0, m * 256, sizex, sizey,
			bandR, sizex, sizey, GDT_Float32, 0, 0);
		poFusDS->GetRasterBand(2)->RasterIO(GF_Write, 0, m * 256, sizex, sizey,
			bandG, sizex, sizey, GDT_Float32, 0, 0); 
		poFusDS->GetRasterBand(3)->RasterIO(GF_Write, 0, m * 256, sizex, sizey,
			bandB, sizex, sizey, GDT_Float32, 0, 0);
	}

	CPLFree(bandR);
	CPLFree(bandG);
	CPLFree(bandB);
	CPLFree(bandH);
	CPLFree(bandS);
	CPLFree(bandP);

	GDALClose(poMulDS);
	GDALClose(poPanDS);
	GDALClose(poFusDS);

	return 0;
}
```

## 四、实验结果

***多光谱图像***

![Mul_large](/img6/Mul_large.PNG)

***全色图像***

![Pan_large](/img6/Pan_large.PNG)

***方式一（256 X 256）结果图像***

![result1](/img6/result1.PNG)

***方式二（图像宽度  X 256）结果图像***

![result2](/img6/result2.PNG)



**第二种方式的处理效率高于第一种方式的处理效率。**



## 五、实验总结

1、每处理一个新的分块时，需要重新设置语句中分块在图像中的起始坐标，以及分块高度和宽度（行或列方向的最后一个分块，其宽度或高度不一定是256）。

2、计算出行或列方向上的分块数，可通过分块数来控制指针的移动。