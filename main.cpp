#include <headers/mpi.h>
#include <iostream>
#include <cstring>
#include <chrono>
#include <math.h>
#include <queue>
#include <fstream>

using namespace std;

//число кратно 10 т.к. числа в десятичной системе, определяет диапозон чисел для перемножения
const int maxNumBase = 10;
//условие, определяющее длину числа, начиная с которой происходит обычное умножение
const int simpleMultNumCond = 1;

//разбивает строку и делает из неё массив чисел, кратных maxNumBase
vector<int> numStrToVec(string str) {
    vector<int> numVec;
    unsigned int digitMultiplier = 1;
    int num = 0;

    //перебираем значения строки
    for(size_t i = 0; i < str.size(); ++i){
        /* формируем число в зависимости от разряда
         * при maxNumBase равном 100, число будет двузначным
         */
        num += (str[i] - '0') * digitMultiplier;
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

//делает размер массива кратным двум степени заданного числа
void fetchVec(vector<int>& vec, size_t size) {
    while (size & (size - 1)) {
        ++size;
    }
    vec.resize(size);
}

//перемножает два массива чисел и возвращает вектор полученных значений
vector<int> multNumVecs(const vector<int> a, const vector<int> b) {
    size_t size = a.size();
    vector<int> ret(size * 2);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            ret[i + j] += a[i] * b[j];
        }
    }

    return ret;
}

//перемножает два массива чисел алгоритмом Карацубы и возвращает вектор полученных значений
vector<int> multNumVecsKaratsuba(const vector<int> A, const vector<int> B) {
    //N - чётное число
    size_t N = A.size();

    //если кол-во чисел в A меньше заданного, выходим из рекурсии, используя обычное умножение
    if (N <= simpleMultNumCond) {
        vector<int> ret = multNumVecs(A, B);
        return ret;
    }

    vector<int> ret(2 * N);
    size_t Ndiv2 = N / 2;

    //число, полученное из старших разрядов x
    vector<int> a (A.begin(), A.begin() + Ndiv2);
    //число, полученное из младших разрядов x
    vector<int> b (A.begin() + Ndiv2, A.end());
    //число, полученное из старших разрядов y
    vector<int> c (B.begin(), B.begin() + Ndiv2);
    //число, полученное из младших разрядов y
    vector<int> d (B.begin() + Ndiv2, B.end());

    ////////вычисляем по формуле A*B=(a+b)(c+d)-a*c-b*d////////

    //перемножаем значения младших разрядов
    vector<int> bd = multNumVecsKaratsuba(b, d);
    //перемножаем значения страших разрядов
    vector<int> ac = multNumVecsKaratsuba(a, c);

    vector<int> aPb(Ndiv2);
    vector<int> cPd(Ndiv2);

    for (size_t i = 0; i < Ndiv2; ++i) {
        aPb[i] = a[i] + b[i];
        cPd[i] = c[i] + d[i];
    }

    vector<int> aPdMcPd = multNumVecsKaratsuba(aPb, cPd);

    for (size_t i = 0; i < N; ++i) {
        aPdMcPd[i] -= ac[i] + bd[i];
    }
    for (size_t i = 0; i < N; ++i) {
        ret[i] = ac[i];
    }
    for (size_t i = N; i < 2 * N; ++i) {
        ret[i] = bd[i - N];
    }
    for (size_t i = Ndiv2; i < N + Ndiv2; ++i) {
        ret[i] += aPdMcPd[i - Ndiv2];
    }

    return ret;
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

//многопоточное умножение, перехватывает управление над подаваемой группой процессов
vector<int> multNumVecsKaratsubaProc(const vector<int> A, const vector<int> B, int lastRank, MPI_Comm lastComm)
{
    //N - чётное число
    size_t N = A.size();

    //ранк текущего процесса
    int currRank;
    MPI_Comm_rank(lastComm, &currRank);

    //если кол-во чисел в A меньше заданного, выходим из рекурсии, используя обычное умножение
    if (N <= simpleMultNumCond) {
        vector<int> ret = multNumVecs(A, B);
        return ret;
    }

    vector<int> ret(N * 2);

    /* проверяем число процессов текущего коммуникатора
     * используя для этого группы
     */
    int currGroupSize;
    MPI_Group currGroup;
    MPI_Comm_group(lastComm, &currGroup);
    MPI_Group_size(currGroup, &currGroupSize);
    MPI_Group_free(&currGroup);

    /* если поток в коммуникторе не один, то распределяем работу
     * иначе используем собственно однопроцессорное перемножение
     */
    if(currGroupSize != 1)
    {
        int newRank;
        int color;

        size_t Ndiv2 = N / 2;

        vector<int> w1, w2;
        vector<int> w3, w4;

        //каждый из коммуникаторов считает свою часть разрядов
        if(lastRank % 2 == 1) {
            w1 = {A.begin(), A.begin() + Ndiv2};
            w2 = {A.begin() + Ndiv2, A.end()};
            w3 = {B.begin(), B.begin() + Ndiv2};
            w4 = {B.begin() + Ndiv2, B.end()};
            color = 1;
        }
        else {
            w1 = {A.begin() + Ndiv2, A.end()};
            w2 = {A.begin(), A.begin() + Ndiv2};
            w3 = {B.begin() + Ndiv2, B.end()};
            w4 = {B.begin(), B.begin() + Ndiv2};
            color = 0;
        }

        ////////вычисляем по формуле A*B=(a+b)(c+d)-a*c-b*d////////
        ////////          a,b,c,d ~ w1,w2,w3,w4            ////////

        vector<int> w1Pw2(Ndiv2);
        vector<int> w3Pw4(Ndiv2);

        //считаем
        for (size_t i = 0; i < Ndiv2; ++i) {
            w1Pw2[i] = w1[i] + w2[i];
            w3Pw4[i] = w3[i] + w4[i];
        }

        //делим коммуникатор на две группы, чтобы процессы получили свою часть работы
        MPI_Comm newComm;
        MPI_Comm_split(lastComm, color, lastRank, &newComm);
        MPI_Comm_rank(newComm, &newRank);

        //входим в рекурисию
        vector<int> Res1_w1Mw3 = multNumVecsKaratsubaProc(w1, w3, newRank, newComm);
        vector<int> w1Pw2Mw3Pw4 = multNumVecsKaratsubaProc(w1Pw2, w3Pw4, lastRank, lastComm);

        //главный(0) поток нового коммуникатора отправляет массив w1Mw3 верхнему в иерархии коммуникатору
        if(lastRank != 0 && newRank == 0) {
            int size = Res1_w1Mw3.size();
            MPI_Send(&size, 1, MPI_INT, 0, 0, lastComm);//размер массива

            MPI_Datatype intArrType;
            MPI_Type_contiguous(size, MPI_INT, &intArrType);
            MPI_Type_commit(&intArrType);

            MPI_Send(&Res1_w1Mw3[0], 1, intArrType, 0, 0, lastComm);//сам массив
            MPI_Type_free(&intArrType);
        }

        //главный(0) поток верхнего в иерархии коммуникатора получает сообщения от подчинённых и заканчивает вычисления
        if(lastRank == 0)
        {
            MPI_Status status;

            int size;

            MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, lastComm, &status);//размер массива
            MPI_Datatype intArrType;
            MPI_Type_contiguous(size, MPI_INT, &intArrType);
            MPI_Type_commit(&intArrType);

            vector<int> Res2_w1Mw3(size);
            MPI_Recv(&Res2_w1Mw3[0], 1, intArrType, MPI_ANY_SOURCE, MPI_ANY_TAG, lastComm, &status);//сам массив

            //считаем...

            for (size_t i = 0; i < N; ++i) {
                w1Pw2Mw3Pw4[i] -= Res2_w1Mw3[i] + Res1_w1Mw3[i];
            }

            //старшие разряды результата
            for (size_t i = 0; i < N; ++i) {
                ret[i] = Res2_w1Mw3[i];
            }
            //младшие разряды результата
            for (size_t i = N; i < 2 * N; ++i) {
                ret[i] = Res1_w1Mw3[i - N];
            }

            for (size_t i = Ndiv2; i < N + Ndiv2; ++i) {
                ret[i] += w1Pw2Mw3Pw4[i - Ndiv2];
            }

            //завершаем работу процесса
            MPI_Comm_free(&newComm);
            MPI_Type_free(&intArrType);

            return ret;
        }

        //завершаем вычисления
        MPI_Comm_free(&newComm);

        return A;
    }
    //если процесс в коммуникаторе один, то это атомарная операция
    else {
        return multNumVecsKaratsuba(A, B);
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
    ifstream file("C:\\Users\\PC\\Documents\\mpi\\nums.txt");
    while(getline(file, str)){
        readBigNums.push_back(str);
    }
    file.close();

    //////////объявляем значения для цикла//////////

    int N = 0;

    string firstNumStr = readBigNums[0];

    string currNumStr = "";

    //формируем первое число
    vector<int> firstNumsVec = numStrToVec(firstNumStr);

    vector<int> secondNumsVec;
    vector<int> res;

    int resSize;

    size_t i = 1;

    //цикл проходит по всем числам
    auto start = std::chrono::steady_clock::now();
    while(i < readBigNums.size()) {
        currNumStr = readBigNums[i];
        secondNumsVec = numStrToVec(currNumStr);

        //приводим числа к одному размеру, кратному 2
        N = max(firstNumsVec.size(), secondNumsVec.size());
        fetchVec(firstNumsVec, N);
        fetchVec(secondNumsVec, N);

        //перемножаем числа в многопроцессорном режиме
        res = multNumVecsKaratsubaProc(firstNumsVec, secondNumsVec, procRank, MPI_COMM_WORLD);
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
        std::cout << endl << "elapsed time " << elapsedTime.count()/1000 << "s\n";

        for (size_t i = 0; i < res.size(); ++i) {
            res[i + 1] += res[i] / maxNumBase;
            res[i] %= maxNumBase;
        }

        printNumsVec(res);
    }
    MPI_Finalize();

    return 0;

}
