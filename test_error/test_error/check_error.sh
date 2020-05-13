list_name=(aims_new  SHEval            sirius_gsl aims_old  SHEval_quadruple  SHEval_single sirius_opt)
for name in ${list_name[*]}
do
python check_error.py error_$name > mae_$name
done
