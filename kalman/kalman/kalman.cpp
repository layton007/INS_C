#include "pch.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include "matrixCalculate.h"
#include "info.h"


using namespace std;



void renewBodyData(double *Cn2b, double V[], double f[], double &latit, double &longti, double &h)
{
	double  Rm, Rn;
	Rm = Re * (1 - 2 * f_ + 3 * f_*sin(latit)*sin(latit)), Rn = Re * (1 + f_ * sin(latit) * sin(latit));

	double tmp1 = V[1] / (Rm + h), tmp2 = 2.0*Wie + V[0] / ((Rn + h)*cos(latit));
	double deta_V[3], fe, fn, fu;

	fe = Cn2b[0] * f[0] + Cn2b[3] * f[1] + Cn2b[6] * f[2];		//载体下比力转东北天
	fn = Cn2b[3] * f[1] + Cn2b[4] * f[1] + Cn2b[7] * f[2];
	fu = Cn2b[6] * f[2] + Cn2b[5] * f[1] + Cn2b[8] * f[2];

	deta_V[0] = (fe + tmp2 * sin(latit)*V[1] - tmp2 * cos(latit) * V[2]) * sampleTime;
	deta_V[1] = (fn - tmp2 * sin(latit)*V[0] - tmp1 * V[2]) * sampleTime;
	deta_V[2] = (fu + tmp2 * cos(latit)*V[0] + tmp1 * V[1] - g_) * sampleTime;

	latit = latit + V[1] * sampleTime / (Rm + h);		//经纬高更新
	longti = longti + V[0] * sampleTime / ((Rn + h) * cos(latit));
	h = h + V[2] * sampleTime;

	for (int i = 0; i < 3; ++i)	//速度更新
		V[i] += deta_V[i];
}

/*
 *Xk 15x1, Fk 15x15, Pk 15x15, Qk 15x15
 *Kk 15x6, Hk  6x15, Rk  6x6,  Zk  6x1
 *Xk为状态矢量，Fk为状态转移矩阵， Pk为状态预测协方差矩阵，Qk表示与INS误差相关的系统噪声的协方差矩阵
 *Kk为卡尔曼增益，Hk为tk时刻的量测矩阵，Rk为量测噪声矩阵， Zk为量测矢量
*/
void kalman(double Xk[], double Fk[], double Pk[], double Qk[],
	double Kk[], double Hk[], double Rk[], double Zk[])
{
	double Xk_1[N] = { 0 }, Pk_1[N*N] = { 0 }, tmp1[N*N] = { 0 }, tmp2[N*N] = { 0 };

	for (int i = 0, j = 0; i < N*N; ++i) {
		if (i % N == 0) {
			Xk_1[j] = Xk[j];
			j++;
		}
		Pk_1[i] = Pk[i];
	}

	//	printMat(Fk, N, N);

		//1 -- 预报状态矢量
	mXm(N, N, 1, Fk, Xk_1, Xk);			//Xk(-) = Fk/k-1 * Xk_1(+)

	//2 -- 预报协方差矩阵
	mXtrm(N, N, N, Pk_1, Fk, tmp1);		//Pk(-) = Fk/k-1 * Pk_1(+) * Fk/k-1` + Qk
	mAddMult(N, N, N, Qk, Fk, tmp1, Pk);

	//3 -- 修正卡尔曼增益
	mXtrm(N, N, M, Pk, Hk, tmp1);		//Kk(-) = Pk(-) * Hk` * {Hk * Pk(-) * Hk` + Rk}^-1;
	mAddMult(M, N, M, Rk, Hk, tmp1, tmp2);

	dcinv(tmp2, M);
	mXm(N, M, M, tmp1, tmp2, Kk);

	//4 -- 修正状态估计
	double tmp_M6x1[M], tmp_Xk[N];					//Xk(-) = Xk(-) + Kk[Zk - Hk * Xk(-)]
	mXm(M, N, 1, Hk, Xk, tmp_M6x1);
	for (int i = 0; i < M; ++i)
		tmp_M6x1[i] = Zk[i] - tmp_M6x1[i];
	mAddMult(N, M, 1, Xk, Kk, tmp_M6x1, tmp_Xk);
	for (int i = 0; i < N; ++i)
		Xk[i] = tmp_Xk[i];

	//5 -- 修正协方差矩阵
	mXm(N, M, N, Kk, Hk, tmp1);		//Pk(+) = [I - Kk * Hk] * Pk(-)
	for (int i = 0; i < N; ++i)
		tmp1[i*N + i] -= 1;
	for (int i = 0; i < N*N; ++i)
		tmp1[i] = -tmp1[i];
	mXm(N, N, N, tmp1, Pk, tmp2);
	for (int i = 0; i < N*N; ++i)
		Pk[i] = tmp2[i];

}

void initFilter(double R0[], double X0[], double P0[], double Hk[])
{
	for (int i = 0; i < 6; ++i)
		Hk[i * 16] = 1;

	double X[] = { 10, 10, 15, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
	for (int i = 0; i < 15; ++i)
		X0[i] = X[i];

	R0[0] = 100;
	R0[7] = 100;
	R0[14] = 225;
	R0[21] = 0.04;
	R0[28] = 0.04;
	R0[35] = 0.04;

	P0[0] = 100;
	P0[16] = 100;
	P0[32] = 225;
	P0[48] = 0.04;
	P0[64] = 0.04;
	P0[80] = 0.04;
	P0[96] = 0.01;
	P0[112] = 0.01;
	P0[128] = 0.01;
	P0[144] = 0.0001;
	P0[160] = 0.0001;
	P0[176] = 0.0001;
	P0[192] = 0.0001;
	P0[208] = 0.0001;
	P0[224] = 0.0001;
}

void renewFk(double F[], double Cn2b[], double f[], double lati, double high, double &Rm, double &Rn)
{
	double sec;
	sec = 1 / cos(lati);
	Rm = Re * (1 - 2 * f_ + 3 * f_*sin(lati)*sin(lati));
	Rn = Re * (1 + f_ * sin(lati)*sin(lati));
	F[4] = 1.0 / (Rm + high);
	F[18] = sec / (Rn + high);
	F[84] = F[4];
	F[88] = -1.0 / (Rn + high);
	F[123] = tan(lati) * F[88];
	F[35] = 1;

	F[52] = f[2], F[53] = -f[1], F[66] = -f[2], F[68] = f[0], F[81] = f[1], F[82] = -f[0];

	F[57] = Cn2b[0], F[58] = Cn2b[3], F[59] = Cn2b[6], F[72] = Cn2b[1], F[73] = Cn2b[4], F[74] = Cn2b[7], F[87] = Cn2b[2], F[88] = Cn2b[5], F[89] = Cn2b[8];
	F[99] = Cn2b[0], F[100] = Cn2b[3], F[101] = Cn2b[6], F[114] = Cn2b[1], F[115] = Cn2b[4], F[116] = Cn2b[7], F[129] = Cn2b[2], F[130] = Cn2b[5], F[131] = Cn2b[8];

	F[144] = -0.005, F[160] = -0.005, F[176] = -0.005, F[192] = -0.005, F[208] = -0.005, F[224] = -0.005;
}

void renewCn2b(double Cn2b[], double p, double r, double y)
{
#if 0
	Cn2b[0][0] = cos(tR)*cos(tH) + sin(tR)*sin(tP)*sin(tH);
	Cn2b[0][1] = -cos(tR)*sin(tH) + sin(tR)*sin(tP)*cos(tH);
	Cn2b[0][2] = -sin(tR)*cos(tP);
	Cn2b[1][0] = cos(tP)*sin(tH);
	Cn2b[1][1] = cos(tP)*cos(tH);
	Cn2b[1][2] = sin(tP);
	Cn2b[2][0] = sin(tR)*cos(tH) - cos(tR)*sin(tP)*sin(tH);
	Cn2b[2][1] = -sin(tR)*sin(tH) - cos(tR)*sin(tP)*cos(tH);
	Cn2b[2][2] = cos(tR)*cos(tP);
#endif

	Cn2b[0] = cos(y) * cos(r) - sin(y) * sin(p) * sin(r);
	Cn2b[1] = -sin(y) * cos(p);
	Cn2b[2] = cos(y) * sin(r) + sin(y) * sin(p) * cos(r);
	Cn2b[3] = sin(y) * cos(r) + sin(p) * cos(y) * sin(r);
	Cn2b[4] = cos(y) * cos(p);
	Cn2b[5] = sin(y) * sin(r) - sin(p) * cos(y) * cos(r);
	Cn2b[6] = -cos(p) * sin(r);
	Cn2b[7] = sin(p);
	Cn2b[8] = cos(p) * cos(r);
}

/*
 *r[0-2] -> latitude longtitude height
 *V[0-2] -> ve vn ve
 *f[0-2] -> fe fn fu
 */
void renewVelocityAndPosition(double r[], double V[], double f[], double Rm, double Rn)
{
	const double sampletime = 1.0;
	double temp1 = V[1] / (Rm + r[2]);	//dlatit
	double temp2 = 2.0*Wie + V[0] / ((Rn + r[2])*cos(r[0]));//2Wie+dlongi
	double Vt[3];

	Vt[0] = V[0] + (f[0] + temp2 * sin(r[0])*V[1] - temp2 * cos(r[0])*V[2])*sampletime;//实际东向速度更新
	Vt[1] = V[1] + (f[1] - temp2 * sin(r[0])*V[0] - temp1 * V[2])*sampletime;//实际北向速度更新
	Vt[2] = V[2] + (f[2] + temp2 * cos(r[0])*V[0] + temp1 * V[1])*sampletime;  //实际天向速度更新

	r[0] += V[1] * sampletime / (Rm + r[2]);//实际的纬度更新
	r[1] += V[0] * sampletime / ((Rn + r[2])*cos(r[0]));//实际的经度更新
	r[2] += V[2] * sampletime;   //this is new added 

	for (int i = 0; i < 3; ++i)
		V[i] = Vt[i];
}

void renewZk(double Zk[], double r_imu[3], double v_imu[3], double r[3], double v[3])
{
	Zk[0] = r_imu[0] - r[0];
	Zk[1] = r_imu[1] - r[1];
	Zk[2] = r_imu[2] - r[2];
	Zk[3] = v_imu[0] - v[0];
	Zk[4] = v_imu[1] - v[1];
	Zk[5] = v_imu[2] - v[2];
}

void writeToFile(ofstream out[], double t_value[])
{
	out[0] << to_string(t_value[0] * H_d) << " ";		//方便与之前的数据对比，得转换成读取之前的格式
	out[1] << to_string(t_value[1] * H_d) << " ";
	out[2] << to_string(t_value[2]) << " ";
	out[3] << to_string(t_value[3] * Km_h) << " ";
	out[4] << to_string(t_value[4] * Km_h) << " ";
	out[5] << to_string(t_value[5] * Km_h) << " ";
	out[6] << to_string(t_value[6] * H_d) << " ";
	out[7] << to_string(t_value[7] * H_d) << " ";
	out[8] << to_string(t_value[8] * H_d) << " ";
	out[9] << to_string(t_value[9]) << " ";
	out[10] << to_string(t_value[10]) << " ";
	out[11] << to_string(t_value[11]) << " ";
	out[12] << to_string(t_value[12] / g_) << " ";
	out[13] << to_string(t_value[13] / g_) << " ";
	out[14] << to_string(t_value[14] / g_) << " ";
}

int main()
{
	double Rm, Rn, r[3], f[3], V[3], Cn2b[9], t_value[15];	//t_value[]分别为纬度、经度、高度、东北天的速度、3轴陀螺仪、3轴加速度
	double &lati = t_value[0], &longti = t_value[1], &height = t_value[2], &ve = t_value[3], &vn = t_value[4]\
		, &vu = t_value[5], &pitch = t_value[6], &roll = t_value[7], &yaw = t_value[8], &wx = t_value[9] = 0, \
		&wy = t_value[10] = 0, &wz = t_value[11] = 0, &ax = t_value[12], &ay = t_value[13], &az = t_value[14];

	double heading, velocity;

	double Xk[N] = { 0 }, Fk[N*N] = { 0 }, Pk[N*N] = { 0 }, Qk[N*N] = { 0 }, Kk[N*M] = { 0 }, Hk[M*N] = { 0 }, Rk[M*M] = { 0 }, Zk[M] = { 0 };
	double rk_1[3], vk_1[3], fk_1[3];	//imu的位置、速度、比力

	ifstream in[11];	//从磁盘读取11个文件
	ofstream out[15];	//

	vector<string> vec_filename = { "./dat/ax.dat" , "./dat/ay.dat" , "./dat/az.dat","./dat/latitude.dat"\
			 ,"./dat/longtitude.dat", "./dat/height.dat","./dat/velocity.dat","./dat/heading.dat"\
			, "./dat/pitch.dat" , "./dat/roll.dat" , "./dat/yaw.dat" };		//读取的文件名

	vector<string> vec_out_filename = { "./out/latitude", "./out/longtitude", "./out/height", "./out/ve", "./out/vn",\
	"./out/vu", "./out/pitch", "./out/roll", "./out/yaw", "./out/wx", "./out/wy", "./out/wz", "./out/ax", "./out/ay", "./out/az" };

	string line;
	vector<double> vec_in[11];	//定义11个顺序容器vector按顺序存放打开的文件
	for (int i = 0; i < 11; ++i) {
		in[i].open(vec_filename[i]);
		line.clear();
		while (getline(in[i], line)) {
			istringstream record(line);
			string word;
			while (record >> word) {
				vec_in[i].push_back(stod(word));
			}
		}
		in[i].close();
	}

	for (int i = 0; i < 15; ++i)
		out[i].open(vec_out_filename[i]);		//创建输出文件


	initFilter(Rk, Xk, Pk, Hk);
	double imu_r[3], imu_v[3] = { 0 };
	imu_r[0] = vec_in[3][0] * D_h, imu_r[1] = vec_in[4][0] * D_h, imu_r[2] = vec_in[5][0];
	for (int cnt = 1, number = vec_in->size(); cnt < number; ++cnt) {	//循环次数由GPS的数据个数决定

		pitch = vec_in[8][cnt] * D_h;
		roll = vec_in[9][cnt] * D_h;
		yaw = vec_in[10][cnt] * D_h;

		lati = vec_in[3][cnt] * D_h;
		longti = vec_in[4][cnt] * D_h;
		height = vec_in[5][cnt];


		ax = (vec_in[0][cnt] - dt_ax) * g_;		//减去零飘
		ay = (vec_in[1][cnt] - dt_ay) * g_;
		az = (vec_in[2][cnt] - dt_az) * g_;
		f[0] = ax, f[1] = ay, f[2] = az;		//当前 a_gps

		heading = vec_in[6][cnt] * D_h;
		velocity = vec_in[7][cnt] * M_s;

		ve = velocity * sin(heading);
		vn = velocity * cos(heading);
		vu = 0.0;

		V[0] = fabs(ve), V[1] = fabs(vn), V[2] = vu;	//当前V_gps

		renewCn2b(Cn2b, pitch, roll, yaw);
		renewFk(Fk, Cn2b, f, lati, height, Rm, Rn);
		renewBodyData(Cn2b, f, imu_v, imu_r[0], imu_r[1], imu_r[2]);

		white(Qk);		//更新Qk
	//	renewVelocityAndPosition(rk_1, vk_1, fk_1, Rm, Rn);
		renewZk(Zk, rk_1, vk_1, r, V);

		kalman(Xk, Fk, Pk, Qk, Kk, Hk, Rk, Zk);


		for (int i = 0; i < N; ++i) {			//每过一次kalman, 修正一次值

		//	cout << Xk[i] << " "; 
			t_value[i] -= Xk[i];
			//	cout << t_value[i] << " ";
		}
		writeToFile(out, t_value);
	}

	for (int i = 0; i < 15; ++i)		//关闭文件
		out[i].close();

	return 0;
}
