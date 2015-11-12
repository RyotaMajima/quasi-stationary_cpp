#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "libfftw3-3.lib")
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<utility>
#include<algorithm>
#include<functional>
#include<complex>
#include<chrono>
#include<fftw3.h>
using namespace std;
using namespace std::chrono;
using Complex = complex<double>;
using vC = vector<Complex>;
using vvC = vector<vC>;
inline fftw_complex* fftwcast(Complex* f){ return reinterpret_cast<fftw_complex*>(f); }

const double X_BEGIN = -5.0, X_END = 25.0; //系の両端
const double L = X_END - X_BEGIN; //空間幅
const int N = 1024; //空間分割数
const double DELTA_X = L / N;

const double T_END = 75; //終了時刻
const int TN = 500; //時間分割数
const double dt = T_END / TN; //時間刻み幅

const double S = 0.9; //波束の幅
const double X = -0.5; //初期波束の原点からのずれ
const double X_OPT = 4.0; //光学ポテンシャルをかける位置
const double b = 1.0 / 3.0; //3次ポテンシャルの係数

const double E_BEGIN = -1.5, E_END = 1.0; //探索するエネルギーの両端
const int EN = 1000; //エネルギー分割数
const double dE = (E_END - E_BEGIN) / EN; //エネルギー刻み幅

double i2x(int i){
    return X_BEGIN + (i + 1) * DELTA_X;
}
double i2k(int i){
    return 2 * M_PI * (i < N / 2 ? i : i - N) / L;
}
double i2E(double begin, int index, double width){
    return begin + index * width;
}

double groundState(double x, double X){
    return pow(1 /(S * S * M_PI), 1.0 / 4.0) * exp(-(x - X) * (x - X) / (2 * S * S));
}
double firstExcited(double x, double X){
    return (sqrt(2) / (S)) * (x - X) * pow(1 / (S * S * M_PI), 1.0 / 4.0) * exp(-(x - X) * (x - X) / (2 * S * S));
}
double secondExcited(double x, double X){
    return sqrt(1.0 / ((3.0 / 2.0) * pow(S, 4) - pow(S, 2) + (1.0 / 2.0))) * (1 / sqrt(2)) * (2 * (x - X) * (x - X) - 1) * pow(1 / (S * S * M_PI), 1.0 / 4.0) * exp(-(x - X) * (x - X) / (2 * S * S));
}

//初期化用関数の定義
void init(vC &f){
    for (int i = 0; i < N; i++){
        f[i] = firstExcited(i2x(i), X);
    }
}

//ポテンシャル関数の定義
double V(double x){
    return (1.0 / 2.0) * x*x - (b / 3.0) * x*x*x - 1.0 / (6 * b * b);
    //return (1.0 / 2.0) * x*x;
}

//光学ポテンシャルの定義
double V_opt(double x){
    return x > X_OPT ? x*x : 0.0;
}

//Simpson積分
double simpson(vC &f){
    double S_even, S_odd;

    S_even = S_odd = 0.0;

    for (int i = 1; i < (N / 2) - 1; i++){
        S_even += norm(f[2 * i]);
    }

    for (int i = 1; i < N / 2; i++){
        S_odd += norm(f[2 * i - 1]);
    }

    return DELTA_X * (norm(f[0]) + 2 * S_even + 4 * S_odd + norm(f[N - 1])) / 3.0;
}

//dtだけ時間発展させる関数
void timeEvolution(vC &f, fftw_plan plan_for, fftw_plan plan_back){
    //ポテンシャル部分の計算
    for (int j = 0; j < N; j++){
        f[j] *= polar(1.0 / N, -V(i2x(j)) * dt); // 1/Nは正規化因子
    }

    //光学ポテンシャル部分の計算
    for (int j = 0; j < N; j++){
        f[j] *= exp(-V_opt(i2x(j)) * dt);
    }

    //f(x, t) -> g(k, t)
    fftw_execute(plan_for);

    //時間発展の計算
    for (int j = 0; j < N; j++){
        f[j] *= polar(1.0, -i2k(j) * i2k(j) * dt / 2.0);
    }

    //f(k, t+dt) -> f(x, t+dt)
    fftw_execute(plan_back);
}

//エネルギーピークのインデックスを求める関数
void getPeaks(vector<pair<double, int>> &peak, vector<double> &res){
    //微分値が正から負に変わったところの値とインデックス
    for (int i = 1; i < EN - 1; i++){
        if (res[i - 1] < res[i] && res[i] > res[i + 1]){
            peak.push_back(make_pair(res[i], i));
        }
    }

    //ピーク値の大きい順にソート
    sort(peak.begin(), peak.end(), [](const pair<double, int> &i, const pair<double, int> &j){ return i.first > j.first; });

    double E_th = 0.02; //しきい値
    //しきい値以下の要素を削除
    peak.erase(remove_if(peak.begin(), peak.end(), [E_th](pair<double, int> pair) {return pair.first < E_th; }), peak.end());

    //得られたピーク値を表示
    cout << "index" << "\t" << "E" << "\t" << "peak value" << endl;
    cout << setprecision(3);
    for (auto pair : peak){
        cout << pair.second << "\t" << i2E(E_BEGIN, pair.second, dE);
        cout << "\t" << pair.first << endl;
    }
}

//固有状態の抽出
void getEigenfunction(vC &phi, vC &f, fftw_plan plan_for, fftw_plan plan_back, double energy){
    for (int i = 0; i <= TN; i++){
        //積分計算
        for (int j = 0; j < N; j++){
            phi[j] += f[j] * polar(dt, energy * (i * dt));
        }

        timeEvolution(f, plan_for, plan_back);
    }

    for (auto &val : phi){
        val /= T_END;
    }
}

int main(){
    auto start = system_clock::now();
    vC f(N);
    vvC A(EN, vC(N));

    //順方向Fourier変換
    fftw_plan plan_for = fftw_plan_dft_1d(N, fftwcast(f.data()), fftwcast(f.data()), FFTW_FORWARD, FFTW_MEASURE);

    //逆方向
    fftw_plan plan_back = fftw_plan_dft_1d(N, fftwcast(f.data()), fftwcast(f.data()), FFTW_BACKWARD, FFTW_MEASURE);

    init(f); //初期条件f(x,0)の設定

    ofstream ofs;
    char filename[50];

    for (int i = 0; i <= TN; i++){
        sprintf(filename, "./output/timeEvo/output%03d.txt", i);
        ofs.open(filename);
        if (!ofs){
            cerr << "file open error!" << endl;
            exit(1);
        }

        for (int j = 0; j < N; j++){
            ofs << i2x(j) << "\t" << norm(f[j]) << "\t" << V(i2x(j)) << endl;
        }

        ofs.close();

        //積分計算
        for (int j = 0; j < EN; j++){
            for (int k = 0; k < N; k++){
                A[j][k] += f[k] * polar(dt, i2E(E_BEGIN, j, dE) * (i * dt));
            }
        }

        //時間発展
        timeEvolution(f, plan_for, plan_back);
    }

    for (auto &vec : A){
        for (auto &val : vec){
            val /= T_END;
        }
    }

    ofs.open("./output/energy.txt");
    if (!ofs){
        cerr << "file open error!" << endl;
        exit(1);
    }

    vector<double> res(EN); //結果格納用配列
    ofs << scientific;
    for (int i = 0; i < EN; i++){
        res[i] = simpson(A[i]);
        ofs << i2E(E_BEGIN, i, dE) << "\t" << res[i] << endl;
    }

    ofs.close();

    vector<pair<double, int>> peak; //ピーク値とインデックスを格納するpair
    getPeaks(peak, res);

    int peakNum = peak.size(); //しきい値以上のピークの数
    vvC phi(peakNum, vC(N));

    for (int i = 0; i < peakNum; i++){
        init(f);
        getEigenfunction(phi[i], f, plan_for, plan_back, i2E(E_BEGIN, peak[i].second, dE));
    }

    vector<vector<double>> ho(3, vector<double>(N));
    for (int i = 0; i < N; i++){
        double x = i2x(i);
        ho[0][i] = norm(groundState(x, 0.0));
        ho[1][i] = norm(firstExcited(x, 0.0));
        ho[2][i] = norm(secondExcited(x, 0.0));
    }

    //再規格化
    for (int i = 0; i < peakNum; i++){
        double sNorm = simpson(phi[i]);
            for (int j = 0; j < N; j++){
                phi[i][j] = norm(phi[i][j]) / sNorm;
            }
    }

    ofs.open("./output/output_phi.txt");
    if (!ofs){
        cerr << "file open error!" << endl;
        exit(1);
    }

    for (int i = 0; i < N; i++){
        ofs << i2x(i) << "\t" << V(i2x(i)) << "\t";
        for (int j = 0; j < peakNum; j++) {
            ofs << real(phi[j][i]) << "\t" << ho[j][i] << "\t";
        }
        ofs << endl;
    }

    ofs.close();

    auto end = system_clock::now();
    auto dur = end - start;
    auto sec = duration_cast<seconds>(dur).count();
    
    cout << endl;
    cout << "execution time : " << sec << "s" << endl;

    return 0;
}