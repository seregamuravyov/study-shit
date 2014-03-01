#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <iostream>
#include <vector>
#include <sys/socket.h>
#include <netdb.h>
#include <poll.h>
#include <sys/wait.h>

const int MAX_BUF = 128;
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
			write(1, buf, 4);
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

int main()
{
    struct addrinfo hints;
    struct addrinfo *servinfo;
    int status;
    memset(&hints, 0, sizeof hints);
    hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE | AI_ALL;
    if ((status = getaddrinfo(NULL, "8008", &hints, &servinfo)) != 0){
        fprintf(stderr, "getaddrinfo error:%s\n", gai_strerror(status));
        exit(1);
    }
    std::cerr << reinterpret_cast<long> (servinfo->ai_next) << std::endl;
    struct addrinfo *servinfo_old = servinfo;
    std::vector<pollfd> pollfds;

    for (int cnt = 0; servinfo != NULL; servinfo = servinfo -> ai_next)
    {
        std::cerr << cnt << std::endl;
        ++cnt;
        int sockfd = socket(servinfo->ai_family,
                servinfo->ai_socktype, servinfo->ai_protocol);
        if (sockfd < 0)
        {
            std::cerr << "sockfd -1" << std::endl;
            continue;
        }
        int yes = 1;
        if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) < 0)
        {
            std::cerr << "setsockopt -1" << std::endl;
            continue;
        }

        if (bind(sockfd, servinfo->ai_addr, servinfo->ai_addrlen) < 0)
        {
            std::cerr << "bind -1" << std::endl;
            perror(" bind error : ");
            continue;
        }
        if (listen(sockfd, 10) < 0)
        {
            std::cerr << "listen -1" << std::endl;
            continue;
        }
        pollfd tmpfd;
        tmpfd.fd = sockfd;
        tmpfd.events = POLLIN | POLLERR;
        tmpfd.revents = 0;
        pollfds.push_back(tmpfd);
    }
    freeaddrinfo(servinfo_old);
    if (pollfds.empty())
    {
        std::cerr << ("pollfds is empty\n") << std::endl;
        exit(1);
    } else
    {
        std::cout << "pollfds.size == " << pollfds.size() << std::endl;
    }
    std::vector<pollfd> clsfds;
    pollfd infd;
    infd.fd = 0;
    infd.events = POLLIN | POLLERR;
    infd.revents = 0;
	clsfds.push_back(infd);
    char buf[MAX_BUF];
	int pos = 0;
	int cnt;
	bool in_dead = false;
	bool out_dead  = false;
    while (1)
    {
        int ready = poll(pollfds.data(), pollfds.size(), 0);
        if (ready) {
			if (ready < 0)
			{
				if (errno == EINTR)
				{
					std::cerr << "poll interupted" << std::endl;
					continue;
				}
				std::cerr << "poll error : " << strerror(errno) << std::endl;
				exit(1);
			}
			for (int i = 0, n = pollfds.size(); i < n; ++i)
			{
				if (pollfds[i].revents & POLLERR)
				{
					pollfds[i].events = 0;
					close(pollfds[i].fd);
					std::cerr << "FAIL: error while listening socket " << std::endl;
					exit(1);
				}
				if (pollfds[i].revents & POLLIN)
				{
					int sockfd = pollfds[i].fd;
					int cfd = accept(sockfd, NULL, NULL);
					printf("accepted, %d\n", cfd);
					if (cfd < 0)
					{
						perror("cfd < 0");
						continue;
					}
					pollfd tmpfd;
					tmpfd.fd = cfd;
					tmpfd.events = POLLOUT | POLLERR;
					tmpfd.revents = 0;
					clsfds.push_back(tmpfd);
				}
			}
        }
        int clsready = poll(clsfds.data(), clsfds.size(), 0);
		if (clsready) {
			if (clsready < 0)
			{
				if (errno == EINTR)
				{
					std::cerr << "clients poll interupted" << std::endl;
					continue;
				}
				std::cerr << "poll error : " << strerror(errno) << std::endl;
				exit(1);
			}
		}
		if (clsfds[0].revents & POLLERR){
			clsfds[0].events = 0;
			close(clsfds[0].fd);
			std::cerr << "FAIL: error while waiting read" << std::endl;
			exit(1);
		}
		if (!(clsfds[0].revents & POLLIN)){
			continue;
		}
		if (pos < MAX_BUF && !in_dead){
			if (cnt = read(clsfds[0].fd, buf + pos, MAX_BUF - pos)){
				if (cnt < 0){
					perror("Error during reading stdin");
					_exit(1);
				}
				pos += cnt;
			} else {
				in_dead = true;
			}
		}
		for (int i = 1, n = clsfds.size(); i < n; ++i){
			if (clsfds[i].revents & POLLERR){
				clsfds[i].events = 0;
				close(clsfds[i].fd);
				std::swap(clsfds[i], clsfds[n-1]);
				i--;
				clsfds.pop_back();
				n--;
			}
			if (clsfds[i].revents & POLLOUT){
				int clfd = clsfds[i].fd;
				if (pos > 0){
					if ((cnt = write(clfd, buf, pos)) > 0){
					} else {
						perror("Error during writing to file");
						exit(1);
					}
				}
			}
		}
		memmove(buf, buf + cnt, pos -= cnt);
	}
}
