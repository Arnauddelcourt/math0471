
export LIB=${HOME}/src/mumps-5.1.2/lib
export INCLUDE=${HOME}/src/mumps-5.1.2/include

if [ -z "${MKLROOT}" ]; then
    :
else
    :
    #LIB=$LIB:${MKLROOT}/lib/intel64
    #INCLUDE=$INCLUDE:${MKLROOT}/include
fi

echo "INCLUDE=${INCLUDE}"
echo "LIB=${LIB}"
