if [[ "$VIRTUAL_ENV" == "/Users/yael/.venv/amuse" ]]
then
     mpirun -n 1 python source/COSMIC_v3.py
else
     echo "please activate amuse"
fi
