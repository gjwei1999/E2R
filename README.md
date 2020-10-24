# E2R

## description

Function EToRConvert in EToRConvert.cc can convert production cuts in energy to cuts in length for electron, gamma, positron, and proton.

## example

example主体来自Geant4 example/basic/B1（任意一个可以运行的Geant4程序都可以）

只需将E2RConvert.hh和相应src/ include/ 放入example/src/   example/include/中

并在example/src/B1RunAction.cc 中插入175行至286行


运行方法：

cd build

cmake ..

make 

./exampleB1 run1.mac
