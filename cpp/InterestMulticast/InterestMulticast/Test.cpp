#include <iostream>
#include <stdlib.h>
#include <time.h>

int main(){
	const int TOTAL_MULTIGROUP=2;
	const int TOTAL_USER = 100;
	srand((unsigned)time(0));

	//pushUM
	int pushUM[TOTAL_MULTIGROUP][TOTAL_USER] = {};
	
	for (int i = 0; i < TOTAL_MULTIGROUP; i++)
	{
		for (int j = 0; j < TOTAL_USER; j++)
		{
			pushUM[i][j] = rand() % 2;
		//	std::cout << "pushUM[" << i << "][" << j << "] is"
		//		<< pushUM[i][j]<<std::endl;
		}
	}

	//interestUM
	double interestUM[TOTAL_MULTIGROUP][TOTAL_USER] = {};

	for (int i = 0; i < TOTAL_MULTIGROUP; i++)
	{
		for (int j = 0; j < TOTAL_USER; j++)
		{
			interestUM[i][j] = (rand() % 1000)/1000.000 ;
		}
	}

	//omegaM
	double omegaM[TOTAL_MULTIGROUP] = {};

	for (int i = 0; i < TOTAL_MULTIGROUP; i++)
	{
		for (int j = 0; j < TOTAL_USER; j++){
			//感兴趣的阈值设置为0.5
			if (pushUM[i][j] == 1 && interestUM[i][j]>=0.5){
				//如果该用户属于当前遍历的多播组iMultiGroup
				//且该用户对多播组iMultiGroup感兴趣
				omegaM[i] = omegaM[i] + 1;
			}
		}
		omegaM[i] = omegaM[i] / 100.00;
		//std::cout << omegaM[i] << std::endl;
	}

	system("pause");
	return 0;
}