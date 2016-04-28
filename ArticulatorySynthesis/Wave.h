#ifndef _WAVE_H_
#define _WAVE_H_



typedef struct{

	char chRIFF[4];            // 'RIFF' (4 Bytes )
	//long nRIFFLen;             // length of the wave file - 8 Bytes (4Bytes)
    int nRIFFLen;

	char chWAVE[4];            // 'WAVE'
	char chFMT[4];             // 'fmt'
	//long nSize;                // 16 Bytes or 18 Bytes
    int nSize;

	short nFormatTag;          // WAVE_FORMAT_PCM
    short nChannels;           // # of channels
    //long  nSamplesPerSec;      // sampling rate
    //long  nAvgBytesPerSec;     // Bytes per second
    int nSamplesPerSec;
    int nAvgBytesPerSec;
    
    short nBlockAlign;         // nChannels * nBits / 8
	short nBits;               // Bits per sample

	char chDATA[4];            // 'data'
	//long nDataLen;             // length of the data
    int nDataLen;
} WAVEFILEHEADER; 


class CWave
{
public:
	CWave(void);
	~CWave(void);

public:
	WAVEFILEHEADER wHeader;
};

#endif