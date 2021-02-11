//--------------------Подключаемые библиотеки--------------------//
#include <headers/mpi.h>
#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char* argv[]) {
    //--------------------Иициализация--------------------//
    int procNum, procRank, recv;
    int m = 2;//количество итераций
    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    //--------------------Главные процесс--------------------//
    if (procRank == 0) {
        printf(" Num of processors: %3d", procNum);
        printf("\n Main process with rank: %3d", procRank);
        int id = 0;

        printf("\n Send message to process: %3d", 1);
        //отправляем сообщение "0" процессу 1
        MPI_Send(&id, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    //--------------------Остальные процессы--------------------//
    else {
        for(int i = 0; i < m; i++) {
            printf("\n Process %3d receiving message... ", procRank);
            printf("\n Current time is: %3d", MPI_Wtime());
            //процесс ождидает сообщение из любого источника
            //первый попавшийся процесс перехватит управление
            MPI_Recv(&recv, 1, MPI_INT, MPI_ANY_SOURCE,
            MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            printf("\n Process recieved massage: %3d", recv);
            printf("\n Current time is: %3d", MPI_Wtime());

            //recv+1 станет отправляемым сообщением
            ++recv;
            //последний процесс в цепочке
            if(procRank == procNum-1) {
                if(m != i+1){
                    //отправляет сообщение recv процессу с рангом на 1 больше своего
                    printf("\n Process sending massage: %3d", recv);
                    printf("\n Current time is: %3d", MPI_Wtime());
                    MPI_Send(&recv, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                }
            }
            //любой другой процесс
            else {
                //отправляет сообщение recv процессу с рангом на 1 больше своего
                printf("\n Process %3d sending massage: %3d", recv);
                printf("\n Current time is: %3d", MPI_Wtime());
                MPI_Send(&recv, 1, MPI_INT, procRank+1, 0, MPI_COMM_WORLD);
            }
        }
    }
    //--------------------Завершение работы MPI--------------------//
    MPI_Finalize();
    return 0;
}
