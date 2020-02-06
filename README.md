# Relative-Mode-Transfer-Method-RMTM-Routine
A routine for RMTM, applied for simple cantilever vabration
#
This routine trys to help understanding the method in the article "Wei H, Li G, Guo P, et al. Effect of Method Type on the Response of Continuum Vibro-Impact[J]. Shock and Vibration, 2019, 2019."
The file "rk4controlyouchuli_jianhua.m" is the main control file.
Copy all the files in the same directory. Then run "rk4controlyouchuli_jianhua.m" in Matlab. You can get a time history plot of displacement of a regular cantilever tip. The RMTM is applied to the regular cantilever, the simplified case is helpful for understanding the approach of the method.
The detailed annotations are in the "rk4controlyouchuli_jianhua.m", one can try to understand the routine referring to the article.
#
the code of R-K Solver in file "rk4sys.m" is from publised book. 
Only for research use. Contact htwei@zzu.edu.cn.
