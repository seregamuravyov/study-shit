#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <ctype.h>
 
const int MAX_BUF = 2048;
volatile bool sigit = false;
 
void transfer(int fin, int fout){
	char buf[MAX_BUF];
	int pos = 0;
	int cnt;
	bool in_dead = false;
	bool out_dead  = false;
 
	while (!sigit && (pos > 0) || !in_dead){
		if (pos < MAX_BUF && !in_dead){
			if (cnt = read(fin, buf + pos, MAX_BUF - pos)){
				if (cnt < 0){
					perror("Error during reading from file");
 
					_exit(1);
				}
				pos += cnt;
 
			} else {
				in_dead = true;
			}
		}
 
		sleep(1);
 
		if (pos > 4){
			char
			write(1, , );
			memmove(buf, buf + 4, pos - 4);
			pos -= 4;
		}
		if (sigit){
			_exit(0);
		}
 
	}
	return;
}
 
void sigint_handler(int sig){
	sigit = true;
}
 
int main(int argc, char **argv) {
 
 	signal(SIGINT, *sigint_handler);
 
 	int filer;
	if ((filer = open("/dev/random", O_RDONLY)) < 0){
		perror("Failed when open random file");
		_exit(1);
	}
 
 
	int fpid1 = fork();
	if (fpid1) {
		close(0);
		int fpid2 = fork();
		if (fpid2 == 0)	{
			transfer(filer, 1);
		}
		waitpid(fpid2, NULL, 0);
	}
	waitpid(fpid1, NULL, 0);
 
 
	return 0;
}