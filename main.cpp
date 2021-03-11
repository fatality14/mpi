//--------------------Подключаемые библиотеки--------------------//
#include <headers/mpi.h>
#include <iostream>
#include <cstring>
#include <chrono>
#include <time.h>

using namespace std;

int main(int argc, char* argv[])
{
    //--------------------Иициализация--------------------//
    int sizeX = 100;

    int procNum, procRank;
    int m = 5;//кол-во итераций

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    //семя для функции rand(), смещено, чтобы было разное у разных процессов
    srand(time(0)+procRank + 314159 * 451 % 42);

    auto start = std::chrono::steady_clock::now();

    //--------------------Инициализация бу--------------------//
    const int numMsg = 2;
    int sendBuff[procNum][numMsg];
    int* recv = new int[numMsg];

    //--------------------Главный цикл--------------------//
    for(int j = 0; j < m; ++j){
        for(int i = 0; i < procNum; ++i) {
            int mes1 = i+1;
            int mes2 = rand() % 100;
            for(int k = 0; k < procNum; ++k){
                sendBuff[k][0] = mes1;
                sendBuff[k][1] = mes2;
            }

            if(procRank == i){
                cout << "processor " << procRank << " send message " << mes2 << endl;
            }

            MPI_Scatter(sendBuff, numMsg, MPI_INT, recv, numMsg, MPI_INT, i, MPI_COMM_WORLD);

            if(recv[0] == procRank) {
                cout << "processor " << procRank << " recv message " << recv[1] << endl;
            }
            if(recv[0] == procNum && procRank == 0) {
                cout << "processor " << procRank << " recv message " << recv[1] << endl;
            }
        }
    }

    delete[] recv;

    //--------------------Завершение работы MPI--------------------//
    MPI_Finalize();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
