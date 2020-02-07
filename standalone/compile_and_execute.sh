#!/bin/bash
g++ -Wall testOpinion.cpp opinion.cpp -lsocp -lcminpack_s_d -lpthread -o testOpinion.x
./testOpinion.x
