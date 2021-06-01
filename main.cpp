#define _USE_MATH_DEFINES

#include "headers/mpi.h"
#include <iostream>
#include <string>
#include <chrono>
#include <math.h>
#include <queue>
#include <fstream>
#include <complex>

using namespace std;

//число кратно 10 т.к. числа в десятичной системе, определяет диапозон чисел для перемножения
const int maxNumBase = 10;

//число PI
const double PI = M_PI;

typedef complex<double> dComplex;
typedef vector<dComplex> vComplex;
typedef vector<int> vInt;
typedef vector<double> vDouble;

//разбивает строку и делает из неё массив чисел, кратных maxNumBase
vInt numStrToVec(string str) {
    vector<int> numVec;
    unsigned int digitMultiplier = 1;
    int num = 0;

    //перебираем значения строки
    for(int c = str.size() - 1; c >= 0; --c) {
        /* формируем число в зависимости от разряда
         * при maxNumBase равном 100, число будет двузначным
         */
        num += (str[c] - '0') * digitMultiplier;
        digitMultiplier *= 10;
        //добавляем число в результат и повторяем для следующего
        if (digitMultiplier == maxNumBase) {
            numVec.push_back(num);
            num = 0;
            digitMultiplier = 1;
        }
    }

    //добавляем последнее число в результат
    if (num != 0) {
        numVec.push_back(num);
    }

    return numVec;
}

//вывод массива значений в виде обычного числа
void printNumsVec(const vector<int>& vec) {
    //итератор вектора
    auto c = vec.crbegin();

    //пропускаем нули
    while (!*c) {
        ++c;
    }

    //выполняется до конца значений вектора
    while (c != vec.crend()) {
        cout << *c++;//вывод
    }
    cout << endl;
}

//перемножает два массива чисел алгоритмом Шёнхаге — Штрассена и возвращает вектор полученных значений
vComplex FourierTransform(vComplex &vec, bool invert) {
    size_t vecSize = vec.size();

    //выход из рекурсии
    if (vecSize == 1) {
        return vec;
    }

    //делим вектор на две части
    vComplex vec0(vecSize/2),  vec1(vecSize/2);
    for (size_t i = 0, j = 0; i < vecSize; i+=2, ++j) {
        vec0[j] = vec[i];
        vec1[j] = vec[i+1];
    }

    //преобразуем первую и вторую части по формуле
    //см. https://en.wikipedia.org/wiki/Discrete_Fourier_transform
    vec0 = FourierTransform (vec0, invert);
    vec1 = FourierTransform (vec1, invert);

    double ang = 2 * PI/vecSize * (invert ? -1 : 1);
    dComplex w(1), wn(cos(ang), sin(ang));

    for (size_t i = 0; i < vecSize/2; ++i) {
        vec[i] = vec0[i] + w * vec1[i];
        vec[i+vecSize/2] = vec0[i] - w * vec1[i];

        if (invert) {
            vec[i] /= 2,  vec[i+vecSize/2] /= 2;
        }

        w *= wn;
    }

    return vec;
}

//многопоточное умножение, перехватывает управление над подаваемой группой процессов
vComplex FourierTransformMPI(vComplex& vec, const vInt& subdims, bool invert, int lastRank, MPI_Comm lastComm) {
    size_t vecSize = vec.size();

    //выходим из рекурсии
    if (vecSize == 1) {
        return vec;
    }

    //число процессв в текущем коммуникаторе
    int procAmount;
    MPI_Comm_size(lastComm, &procAmount);

    /* если процесс в коммуникторе не один, то распределяем работу
     * иначе используем собственно однопроцессорное перемножение
     */
    if(procAmount > 1) {
        //преобразуем первую и вторую части по формуле
        //см. https://en.wikipedia.org/wiki/Discrete_Fourier_transform

        //первая часть
        vComplex vec0 (vecSize/2);

        if(lastRank < procAmount / 2) {
            for (size_t i = 0, j = 0; i < vecSize; i+=2, ++j) {
                vec0[j] = vec[i];
            }
        }
        else {
            for (size_t i = 1, j = 0; i <= vecSize; i+=2, ++j) {
                vec0[j] = vec[i];
            }
        }

        //делим решётку пополам
        MPI_Comm newComm;
        MPI_Cart_sub(lastComm, &subdims[0], &newComm);

        //ранк процесса в новом коммуникаторе
        int newRank;
        MPI_Comm_rank(newComm, &newRank);

        //входим в рекурсию, преобразуем свою часть вектора
        vec0 = FourierTransformMPI(vec0, subdims, invert, newRank, newComm);

        //главный(0) поток нового коммуникатора отправляет вектор vec0 верхнему в иерархии коммуникатору
        if (lastRank != 0 && newRank == 0) {
            int size = vec0.size();

            vDouble sendVec(size * 2);
            for(int i = 0; i < size; i++) {
                sendVec[i] = vec0[i].real();
                sendVec[i+size] = vec0[i].imag();
            }

            MPI_Send(&size, 1, MPI_INT, 0, 0, lastComm);//размер вектора

            MPI_Datatype longArrType;
            MPI_Type_contiguous(size * 2, MPI_DOUBLE, &longArrType);
            MPI_Type_commit(&longArrType);

            MPI_Send(&sendVec[0], 1, longArrType, 0, 0, lastComm);//сам вектор
            MPI_Type_free(&longArrType);
        }

        //главный(0) поток верхнего в иерархии коммуникатора получает сообщения от подчинённых и заканчивает вычисления
        if(lastRank == 0) {
            MPI_Status status;

            int size;

            MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, lastComm, &status);//размер массива

            MPI_Datatype longArrType;
            MPI_Type_contiguous(size * 2, MPI_DOUBLE, &longArrType);
            MPI_Type_commit(&longArrType);

            vDouble recvVec(size * 2);
            MPI_Recv(&recvVec[0], 1, longArrType, MPI_ANY_SOURCE, MPI_ANY_TAG, lastComm, &status);//сам массив

            MPI_Type_free(&longArrType);

            //заполняем вектор комплексных чисел, полученными значениями
            vComplex vec1(size);
            for(int i = 0; i < size; i++) {
                dComplex v(recvVec[i],recvVec[i+size]);
                vec1[i]=v;
            }

            //считаем...

            double ang = 2 * PI/vecSize * (invert ? -1 : 1);

            dComplex w(1), wn(cos(ang), sin(ang));

            for (size_t i = 0; i < vecSize/2; ++i) {
                vec[i] = vec0[i] + w * vec1[i];
                vec[i+vecSize/2] = vec0[i] - w * vec1[i];

                if (invert) {
                    vec[i] /= 2,  vec[i+vecSize/2] /= 2;
                }

                w *= wn;
            }
        }

        //завершаем вычисления
        MPI_Comm_free(&newComm);
        return vec;
    }
    //если процесс в коммуникаторе один, то это атомарная операция
    else {
        return FourierTransform (vec,invert);
    }
}

void multNumVecsSchonhageStrassenMPI (const vInt & A, const vInt & B, vector<int> & res, const vector<int> &subdims, int lastRank, MPI_Comm lastComm) {
    //переводим числа в комплексные
    vComplex a (A.begin(), A.end());
    vComplex b (B.begin(), B.end());

    //ищем размер результата как 2^n, чтобы точно хватило на все значения
    size_t resSize = 1;
    while (resSize < max (A.size(), B.size())) {
        resSize *= 2;
    }
    resSize *= 2;

    a.resize(resSize);
    b.resize(resSize);

    //вычисляем в многопроцессорном режиме многочлены 'a' и 'b' и заносим его в эти же переменные - прямое преобразование Фурье
    FourierTransformMPI(a, subdims, false, lastRank, lastComm);
    FourierTransformMPI(b, subdims, false, lastRank, lastComm);

    //магия перемножения многочленов
    for (size_t i=0; i<resSize; ++i) {
        a[i] *= b[i];
    }

    //находим получившееся число - обратное преобразование Фурье
    FourierTransform(a, true);

    res.resize(resSize);

    //переводим число из комплексного в обычное
    for (size_t i = 0; i < resSize; ++i) {
        res[i] = int(a[i].real() + 0.5);
    }

    //определяем выход за пределы массива, если такие были
    int rest = 0;
    for (size_t i = 0; i < resSize; ++i) {
        res[i] += rest;
        rest = res[i] / 10;
        res[i] %= 10;
    }
}

int main(int argc, char* argv[])
{
    int procNum, procRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    //считываем числа из файла построчно
    vector<string> readBigNums;
    string str;

    //НУЖНО СОЗДАТЬ ФАЙЛ И ДОБАВИТЬ В НЕГО ЧИСЛА, ОДНА СТРОКА - ОДНО ЧИСЛО, ПУСТЫЕ СТРОКИ НЕЛЬЗЯ!
    ifstream file("C:\\Users\\1234\\Documents\\prjcts\\mpii\\nums.txt");

    while(getline(file, str)){
        readBigNums.push_back(str);
    }
    file.close();

    if((double)((int)log2(procNum)) == log2(procNum)) {
        //////////объявляем значения для цикла//////////

        int resolution = (int)log2(procNum);

        string firstNumStr = readBigNums[0];
        string currNumStr = "";

        //формируем первое число
        vector<int> firstNumsVec = numStrToVec(firstNumStr);
        vector<int> secondNumsVec;
        vector<int> res;

        int resSize;

        size_t i = 1;

        //////////объявляем коммуникатор решётки 2*2*...*2*1//////////

        MPI_Comm gridComm;

        int reorder = 1;
        vInt dims(resolution + 1, 2);
        dims[resolution]=1;
        vInt periods(resolution + 1, 1);
        vInt subdims(resolution+1,1);
        subdims[0]=0;

        MPI_Cart_create(MPI_COMM_WORLD, resolution + 1, &dims[0], &periods[0], reorder, &gridComm);

        //цикл проходит по всем числам
        auto start = std::chrono::steady_clock::now();
        while (i < readBigNums.size()) {
            currNumStr =readBigNums[i];
            secondNumsVec = numStrToVec(currNumStr);

            //перемножаем числа
            multNumVecsSchonhageStrassenMPI(firstNumsVec, secondNumsVec, res, subdims, procRank, gridComm);

            resSize = res.size();

            //после перемножения сообщаем главному процессу результат
            MPI_Bcast(&resSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(procRank != 0) {
                res.resize(resSize);
            }
            MPI_Bcast(&res[0], resSize, MPI_INT, 0, MPI_COMM_WORLD);
            firstNumsVec = res;

            ++i;
        }
        auto end = std::chrono::steady_clock::now();

        //главный процесс правит полученный результат в соответствии с диапозоном значений maxNumBase и выводит его
        std::chrono::duration<double> elapsedTime = end - start;
        if(procRank == 0) {
            std::cout << "elapsed time " << elapsedTime.count()/1000 << "s\n";

            printNumsVec(res);
        }

        MPI_Finalize();
    }

    return 0;

}
