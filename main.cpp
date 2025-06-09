#ifdef ONLINE_JUDGE

#include <cstring>
#include <cstdio>
#include "Point.h"
#include "Strategy.h"

void send(const char *msg) {
    size_t len = strlen(msg);
    printf("%c%c%c%c%s", (unsigned char)(len >> 24), (unsigned char)(len >> 16), (unsigned char)(len >> 8), (unsigned char)(len), msg);
    fflush(stdout);
}

int main()
{
    int noX, noY;
    int lastX, lastY;
    int row, col;
    int* board; //chess board
    int* top; //available position

    scanf("%d %d %d %d", &row, &col, &noX, &noY);
    board = new int[row * col];
    top = new int[col];
    while(true)
    {
        scanf("%d %d", &lastX, &lastY);
        for(int i = 0; i < col; i ++) scanf("%d", &top[i]);
        for(int i = 0; i < row * col; i ++) scanf("%d", &board[i]);

        Point* point = getPoint(row, col, top, board, lastX, lastY, noX, noY);
        char msg[6];
        sprintf(msg, "%d %d", point->x, point->y);
        send(msg);
        clearPoint(point);
    }
    return 0;
}

#endif  // ONLINE_JUDGE
