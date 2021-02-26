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

    int allSum = 0;

    auto start = std::chrono::steady_clock::now();

    //--------------------Запуск цикла--------------------//
    if (procRank == 0) {
        recv = 0;

        cout << "Num of processors: " << procNum << endl;
        cout << "Main process with rank: " << procRank << endl;
        cout << "Start broadcasting" << endl;
        MPI_Bcast(&recv, 1, MPI_INT, procRank, MPI_COMM_WORLD);
    }
    //--------------------Главный цикл--------------------//
    while(true) {
        int globalSum = 0;
        int localSum = 0;
        for(int i = 0; i < sizeX; ++i){
            localSum += rand() % 100;
        }

        cout << "Local sum of process " << procRank << " is " << localSum << endl;
        MPI_Reduce(&localSum, &globalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if(procRank == 0){
            allSum += globalSum;
            cout << "Global sum is " << globalSum << endl;
        }

        //Вычисляем значение root, из которого будем получать сообщение
        int fromProc = (procNum+procRank-1) % procNum;
        cout << "Processor " << procRank << " recv the message from root " <<  fromProc << endl;

        //Получаем сообщение
        MPI_Bcast(&recv, 1, MPI_INT, fromProc, MPI_COMM_WORLD);
        cout << "Processor " << procRank << " recv the message " << recv << endl;

        //Увеличиваем значение счётчика на 1
        int nextCount = recv + 1;
        cout << "Processor " << procRank << " send the message " << nextCount << endl;

        //Потом отправляем сообщение
        MPI_Bcast(&nextCount, 1, MPI_INT, procRank, MPI_COMM_WORLD);

        //Проверяем завершение цикла, m кол-во итераций, а (recv-1)/procNum текующая итерация
        if((recv-1)/procNum == m-1){
            break;
        }
    }

    if(procRank == 0){
        cout << "Final sum is " << allSum << endl;
    }

    //--------------------Завершение работы MPI--------------------//
    MPI_Finalize();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
