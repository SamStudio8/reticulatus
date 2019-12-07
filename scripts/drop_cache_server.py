import socket
import subprocess

HOST = '127.0.0.1'
PORT = 34837
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((HOST, PORT))
s.listen()

while True: 
    conn, addr = s.accept()
    with conn:
        print('Connected by', addr)
        while True:
            data = conn.recv(1024)
            if not data:
                break
            conn.sendall(data)
            print(subprocess.run("sync", shell=True, check=True))
            print(subprocess.run("echo 3 > /proc/sys/vm/drop_caches", shell=True, check=True))

