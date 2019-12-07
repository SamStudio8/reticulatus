import socket

HOST = '127.0.0.1'
PORT = 34837

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((HOST, PORT))
s.sendall(b'DROP')
data = s.recv(1024)
if data:
    print('[snakeshell] Successfully dropped system cache')
