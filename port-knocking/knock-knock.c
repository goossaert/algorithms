/******************************************************
 *  knock-knock 0.1 - GPL License v3.0
 *  Copyright (c) 2010 Emmanuel Goossaert
 *  emmanuel at goossaert dot com
 *
 *  Knock-knock is a daemon that implements port knocking
 *  in order to protect the access to a specific service
 *  on a server. It adds an ip address in the list of the
 *  authorized addresses for a service each time this ip
 *  address knocks on a specified port. After a certain
 *  time, this address is removed from the list: this is
 *  a temporary exception.
 *
 *  The main purpose is to allow just one ip to access
 *  a service during a short period of time, reducing
 *  the risks of brute force attacks. A scan detection
 *  has been included, and prevents the access to be
 *  open if someone knocks sequentially on all ports
 *  during a scan.
 *
 ******************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <syslog.h>
#include <signal.h>

#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>


#define SIZE_LINE        512
#define SIZE_BUFFER      8192


// Parameters have been hardcoded for security reasons. Indeed, a simple 'ps'
// can reveal the parameters passed to a program, and would therefore reveal
// the knocking ports.

// name of the service to be protected
const char *       SERVICE         = "sshd";                  

// port number of the first knock
const unsigned int PORT_KNOCK1     = 8000;                    

// port number of the second knock (can be any value except
// PORT_KNOCK1 - 1 and PORT_KNOCK1 + 1, due to scan detection)
const unsigned int PORT_KNOCK2     = 9000;                    

// timeout delay between the two knocks, in seconds
const unsigned int TIMEOUT_KNOCK   = 5;                       

// timeout delay for the connection after the two knocks, in seconds
const unsigned int TIMEOUT_CONNECT = 5;                       

// path to the network permission file
const char *       FILE_RULE       = "/etc/hosts.allow";      

// path to a temporary file
const char *       FILE_RULETEMP   = "/tmp/knock-knock";      

// name used by the daemon to log its events
const char *       DAEMON_NAME     = "knock-knock";           


struct tcpheader
{
    unsigned short int th_sport;
    unsigned short int th_dport;
    unsigned int th_seq;
    unsigned int th_ack;
    unsigned char th_x2:4, th_off:4;
    unsigned char th_flags;
    unsigned short int th_win;
    unsigned short int th_sum;
    unsigned short int th_urp;
}; // total tcp header length: 20 bytes (=160 bits)


struct ipheader
{
    unsigned char ip_hl:4, ip_v:4;
    unsigned char ip_tos;
    unsigned short int ip_len;
    unsigned short int ip_id;
    unsigned short int ip_off;
    unsigned char ip_ttl;
    unsigned char ip_p;
    unsigned short int ip_sum;
    unsigned int ip_src;
    unsigned int ip_dst;
}; // total ip header length: 20 bytes (=160 bits) 


char* make_rule_line( char * ip )
{
    static char line[ SIZE_LINE ];
    sprintf( line, "%s: %s : allow\n", SERVICE, ip );
    return line;
}


void add_ip( char * ip )
{
    char * line = make_rule_line( ip );
    FILE * file = fopen( FILE_RULE, "at");
    fwrite( line, 1, strlen( line ), file );
    fclose( file ); 
}


void remove_ip( char * ip )
{
    char line[ SIZE_LINE ];
    FILE *file_old = fopen( FILE_RULE, "rt");
    FILE *file_new = fopen( FILE_RULETEMP, "wt");
    
    while( fgets ( line, SIZE_LINE, file_old ) != NULL )
    {
        if( strcmp( make_rule_line( ip ), line ) != 0 )
            fwrite( line, 1, strlen( line ), file_new );
    }
    
    fclose( file_new ); 
    fclose( file_old );
  
    rename( FILE_RULETEMP, FILE_RULE ); 
}



void signal_handler( int sig )
{
 
    switch(sig)
    {
        case SIGHUP:
            syslog(LOG_WARNING, "Received SIGHUP signal.");
            exit(EXIT_FAILURE);
            break;
        
        case SIGTERM:
            syslog(LOG_WARNING, "Received SIGTERM signal.");
            exit(EXIT_FAILURE);
            break;
        
        default:
            syslog(LOG_WARNING, "Unhandled signal (%d) %s", sig, strsignal(sig) );
            break;
    }
}



int main(void)
{
    int fd = socket(AF_INET, SOCK_RAW, IPPROTO_TCP);
    char buffer[SIZE_BUFFER];
    struct ipheader *ip = (struct ipheader *) buffer;
    struct tcpheader *tcp = (struct tcpheader *) (buffer + sizeof(struct ipheader));

    pid_t pid, sid;

    // Error with PORT_KNOCK2 value
    if( PORT_KNOCK2 == PORT_KNOCK1 + 1 || PORT_KNOCK2 == PORT_KNOCK1 - 1 ) {
        exit(EXIT_FAILURE);
    }
 
    // Fork off the parent process
    pid = fork();
    if (pid < 0) {
        exit(EXIT_FAILURE);
    }
    
    // If we got a good PID, then we can exit the parent process.
    if (pid > 0) {
        exit(EXIT_SUCCESS);
    }
 
    // Create a new SID for the child process
    sid = setsid();
    if (sid < 0) {
        // Log the failure
        exit(EXIT_FAILURE);
    }

    setlogmask( LOG_UPTO(LOG_INFO) );
    openlog( DAEMON_NAME, LOG_CONS, LOG_USER );

    syslog( LOG_INFO, "daemon starting up" );  

    signal(SIGHUP, signal_handler);
    signal(SIGTERM, signal_handler);
    signal(SIGINT, signal_handler);
    signal(SIGQUIT, signal_handler);


    // Listen to the network and wait for someone to knock!
    while( read(fd, buffer, SIZE_BUFFER) > 0 )
    {
        char ipaddr[ 16 ];
        
        if( htons(tcp->th_dport) == PORT_KNOCK1 )
        {
            time_t timestamp = time( NULL );
            struct in_addr addr;
            addr.s_addr = ip->ip_src; 

            // Waiting for the confirmation packet
            while( read(fd, buffer, SIZE_BUFFER) > 0 )
            {
                // time is over!
                if( time( NULL ) > timestamp + TIMEOUT_CONNECT )
                    break;

                // this packet is not from the same source 
                if( addr.s_addr != ip->ip_src )
                    continue; 
                
                // detect scanning 
                if( htons(tcp->th_dport) == PORT_KNOCK1 + 1 || htons(tcp->th_dport) == PORT_KNOCK1 - 1 )
                    break;

                // we have an authorized access
                if( htons(tcp->th_dport) == PORT_KNOCK2 )
                {
                    strcpy( ipaddr, inet_ntoa( addr ) );
                    add_ip( ipaddr );
                    syslog( LOG_INFO, "giving access to: %s", ipaddr );  
                    sleep( TIMEOUT_KNOCK ); 
                    remove_ip( ipaddr ); 
                    break;
                }
            }
        }
    }
    
    syslog( LOG_INFO, "exiting, good bye!" );  
    return 0;
}
