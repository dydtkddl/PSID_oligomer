# LAMMPS 소스 디렉토리로 이동
cd ~/downloads/lammps

# ReaxFF 패키지 포함하여 재컴파일
mkdir -p build-reaxff
cd build-reaxff
cmake -D PKG_REAXFF=yes -D CMAKE_BUILD_TYPE=Release ../cmake
make -j 4

# 새로운 실행파일 경로 업데이트
export PATH=~/downloads/lammps/build-reaxff:$PATH
