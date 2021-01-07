#!/bin/bash

> timings_out.txt

echo 44100
echo 44100 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav  44100 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 48000
echo 48000 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav  48000 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 88200
echo 88200 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav  88200 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 96000
echo 96000 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav  96000 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 176400
echo 176400 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav 176400 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 192000
echo 192000 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav 192000 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 352800
echo 352800 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav 352800 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

echo 384000
echo 384000 >> timings_out.txt
(time ./gebouw/resampler data/PinkFloyd-EmptySpaces.wav 384000 &>> timings_out.txt) &>> timings_out.txt
echo "" >> timings_out.txt
echo "" >> timings_out.txt

