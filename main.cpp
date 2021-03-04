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

    int procNum, procRank, recv;
    int m=5;//кол-во итераций

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    //семя для функции rand(), смещено, чтобы было разное у разных процессов
    srand(time(0)+procRank + 314159 * 451 % 42);

    auto start = std::chrono::steady_clock::now();

    //--------------------Запуск цикла--------------------//
    if (procRank == 0) {
        cout << "Num of processors: " << procNum << endl;
    }
    //--------------------Главный цикл--------------------//
    for(int j = 0; j < m; ++j){
        for(int i = 0; i < procNum; ++i) {
            int globalSum = 0;
            int localSum = 0;

            localSum = rand()%100;

            //cout << "Local sum of process " << procRank << " is " << localSum << endl;
            MPI_Reduce(&localSum, &globalSum, 1, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);

            if(globalSum != 0){
                cout << "Process " << procRank << " recieved " << globalSum << endl;
            }
        }
    }

    //--------------------Завершение работы MPI--------------------//
    MPI_Finalize();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
