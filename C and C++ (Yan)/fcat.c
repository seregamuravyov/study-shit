#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>

#define BUF 128

int curf = 0,
    n, wf;
char* wrf, *buf;
int* rf;

void corr_exit(int m) {
    free(buf);
    for (int i = 0; i < m; ++i) {
        close(rf[i]);
    }
    close(wf);
    free(rf);
    exit(0);
}

void incorr_exit(int m) {
    unlink(wrf);
    corr_exit(m);
}

void sig_handler(int sig_no) {
    if (sig_no == SIGUSR1) {
        incorr_exit(n);
    } else if (sig_no = SIGINT) {
        ++curf;
    }
}

int main(int argc, char* arg[]) {
    if (argc == 0) {
        printf("usage: fcat [files_from ..] file_to\n");
        return 0;
    }
 
    signal(SIGINT, sig_handler);
    signal(SIGUSR1, sig_handler); 
    wf = open(arg[argc - 1], O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
    if (wf == -1) {
        printf("can't open target file\n");
        exit(0);
    } 
    n = argc - 2;
    wrf = arg[argc - 1];
    rf = malloc((argc - 2) * sizeof(int));
    for (int i = 1; i < argc - 1; ++i) {
        rf[i - 1] = open(arg[i], O_RDONLY);
        if (rf[i - 1] == -1) incorr_exit(i - 1);
    }
  
    buf = malloc(BUF * sizeof(char));
    int rb, tmp;

    while (curf < argc - 2) {
        if (rb = read(rf[curf], buf, BUF)) {
            if (rb == -1) incorr_exit(n);
            int w = 0;
            while (w < rb) {
                tmp = write(wf, buf + w, rb - w);
                if (tmp == -1) incorr_exit(n);
                w += tmp;
            }
        } else {
            close(curf);
            curf++;
        }
    }  
    corr_exit(n);
}
