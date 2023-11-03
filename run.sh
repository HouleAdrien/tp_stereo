g++ -c main.cpp -o main.o
g++ -o main main.o -lpthread -lX11
./main TurtleD.tif TurtleG.tif