// ArticulatorySynthesis.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"

//
//  main.cpp
//  ArticulatorySynthesis
//
//  Created by fangqiang on 11-10-23.
//  Copyright 2011�� __MyCompanyName__. All rights reserved.



#include "VT_sim.h"
#include "PhysicalConstants.h"

#include "DSP.h"
#include <iostream>
#include "Wave.h"

using namespace std;


int main(int argc, const char * argv[])
{
	//double fs = 1.6e4;
	//CVT_sim vt_s( fs ) ;
	//char *fn = "E:\\a.txt";

	//vt_s.loadconfig( fn, VT );
	//
	//for( int i = 0; i < vt_s.nvt; i++ ){
	//	printf("%d %f %f\n", i, vt_s.vt[i].a[0],  vt_s.vt[i].x[0]);
	//}

	//double a[5];
	//Complex r[4];

	//r[0].re = r[1].re = r[2].re = r[3].re = 1.0;
	//r[0].im = 2; r[1].im = -2;
	//r[2].im = 4; r[3].im = -4;

	//// poly( r, a, 4 );

	double fs = 4.4e4;
	DSP dsp(fs);
	//double wc[1] = {1.0e3};
	//int N = 3;

	//double *a = new double[N+1];
	//double *b = new double[N+1];

	//dsp.butterworth( N, wc, 1, fs, 1, a, b ); 

	//for( int i = 0; i <= N; i++ ){
	//	printf("a[%d]: %f  b[%d]:%f\n", i, a[i], i, b[i] )  ;
	//}

	//delete [] a;
	//delete [] b;

	CVT_sim vt_sim(fs);
	CWave wave;
	double a[100];   // area function of vocal tract �������������
	double l[100];   // length of short tubes �ܵ��ĳ���
	double Ag = 0.0; // cross-sectional area of glottis ���ŵĺ�����

	int i = 0, j = 0; // k = 0;

	FILE *fp;
	char fn_wav[100] = "./wav/s.wav";
	char fn[100] = "./config/Area_s.txt";

	// load the initial shape of vocal tract ���������ĳ�ʼ����״
	fp = fopen(fn, "r");
	if (fp == NULL) {
		cout << "Can not open file: " << fn << endl;
		return 0;
	}

	i = 0;
	while (fscanf(fp, "%lf", &a[i]) != EOF) {
		fscanf(fp, "%lf", &l[i]);
		a[i] = a[i] / 10000;
		l[i] = l[i] / 100;
		i = i + 1;
	}
	fclose(fp);

	/*
	printf("Size of char: %lu\n", sizeof( char) );
	printf("Size of short: %lu\n", sizeof( short ));
	printf("Size of long int: %lu\n", sizeof( long int ) );
	printf("Size of int:  %lu\n", sizeof( int ) );
	printf("Size of word: %lu\n", sizeof( WORD_BIT) );
	*/

	vt_sim.nPos_th =19.25;//���������ݵ�λ��
	vt_sim.Init(a, l, i);//��ʼ�����ùܵ� iΪ�ܵĳ���

	// initialization
	vt_sim.ngl = 1;//���ų���
	vt_sim.InitResidue(&Ag, 1, 8e-2*H2O);//��ʼ���в�

	// set the wave file header ����wave�ļ�ͷ��ʽ
	strcpy(wave.wHeader.chRIFF, "RIFF");//ÿ��wave�ļ���ͷ�ĸ��ֽڱ��ǡ�RIFF��
	wave.wHeader.nRIFFLen = 0;

	strcpy(wave.wHeader.chWAVE, "WAVE");
	strcpy(wave.wHeader.chFMT, "fmt ");
	wave.wHeader.nSize = 16;

	wave.wHeader.nFormatTag = 1; //WAVE_FORMAT_PCM;
	wave.wHeader.nChannels = 1;
	wave.wHeader.nSamplesPerSec = (int)fs;//ȡ��Ƶ��
	wave.wHeader.nBits = sizeof(short) * 8;
	wave.wHeader.nAvgBytesPerSec = wave.wHeader.nBits * wave.wHeader.nChannels * wave.wHeader.nSamplesPerSec / 8;//ÿ���ֽ�
	wave.wHeader.nBlockAlign = wave.wHeader.nChannels * wave.wHeader.nBits / 8;

	strcpy(wave.wHeader.chDATA, "data");
	wave.wHeader.nDataLen = 0;//��ʼ��data����Ϊ0


	// calculate the glottal area�����������
	double a_nv[3] = { 0.0, 0.0, 0.0 };
	double T0 = 1.0 / 90;
	//double t  = 0.0;
	double OQ = 0.75;
	double SQ = 2.0;
	double Ap = 1.0e-4;
	double Ps = 8.0e-2 * H2O;


	int nsteps = (int)(T0*fs);

	long nDataLen = 0;
	double *sp_f = new double[nsteps];
	short *sp_i = new short[nsteps];

	// int info = 0;


	if (!(fp = fopen(fn_wav, "wb"))) {
		cout << "Cann't create wave file: " << fn_wav << endl;
		return 0;
	}

	fwrite(&wave.wHeader, sizeof(WAVEFILEHEADER), 1, fp);//д�������wav�ļ�ͷ

	for (i = 0; i < 20; i++) {
		//memset( sp_f, 0.0, nsteps*sizeof(double));
		//memset( sp_i, 0, nsteps*sizeof(short) );
		for (j = 0; j < nsteps; j++) {
			Ag = vt_sim.glottalArea_dy(T0, j / fs, OQ, SQ, Ap);//������������Ķ�̬��� ����Ag_dy     
			//			cout<<"Ag: "<<Ag<<endl;
			sp_f[j] = vt_sim.td_sim(&Ag, a, l, a_nv, Ps);//����ϵͳ��ʱ��ģ�� ���ø����ܵ�����ѧԪ�� acoust_elm_t()  �޸�����modifyResistance() ���㵥��������Դ ż��������Դ   
														//��������ϵͳ�����Ե�ʽ ����ÿ���ܵ��Ŀ����� �촽ĩβPl, �ǿ�Pn ������Pw���ķ��� �������ܵ��ļ�����Ϣ
			sp_i[j] = (short)(sp_f[j] * 20000 / 2.0);

			//sp = sp_i[j];
			// fwrite( &sp, sizeof(short), 1, fp );
			//fwrite( &sp_i[j], sizeof(short), 1, fp );
		}
		fwrite(sp_i, sizeof(short), nsteps, fp);//д��
		//fwrite( sp_f, sizeof(float), nsteps, fp );
		nDataLen = nDataLen + nsteps * wave.wHeader.nBlockAlign;
	}
	//fclose( fp );


	wave.wHeader.nDataLen = (int)nDataLen;

	wave.wHeader.nRIFFLen = (int)nDataLen + 44 - 8;

	rewind(fp);//ʹ�ļ�fp��λ��ָ��ָ���ļ���ʼ
	fwrite(&wave.wHeader, sizeof(WAVEFILEHEADER), 1, fp);//����д���ļ�ͷ

	fclose(fp);

	delete[] sp_f;
	delete[] sp_i;

	return 0;
}




