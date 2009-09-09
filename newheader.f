      subroutine bodevect(%%vectinputs%%)
%%vectdefs%%
      double complex bode
      DO i=1,n0
         outvect(i) = bode(%%vectcallingscalar%%)
      ENDDO
      END

      subroutine invbodevect(%%vectinputs%%)
%%vectdefs%%
      double complex invbode
      DO i=1,n0
         outvect(i) = invbode(%%vectcallingscalar%%)
      ENDDO
      END

      double complex function invbode(%%scalarinputs%%)
%%scalardefs%%
Cf2py intent(out) invbode
      double complex bode
      invbode = 1.0/bode(%%scalarcallingscalar%%)
      RETURN
      END

      double complex function zcosh(z)
      double complex z
      zcosh = 0.5*(exp(z)+exp(-z))
      RETURN
      END
      
      double complex function zsinh(z)
      double complex z
      zsinh = 0.5*(exp(z)-exp(-z))
      RETURN
      END

      double complex function poly(s,coeffs,n)
Cf2py intent(in) s, coeffs
Cf2py intent(out) poly
Cf2py integer intent(hide),depend(coeffs) :: n = len(coeffs)
      double complex s
      integer n
      double precision coeffs(n)
      poly=0
      DO i=1,n
!         PRINT *,'i=', i
         poly = poly+coeffs(i)*s**(n-i)
      ENDDO
      RETURN
      END

      double complex function bode(%%scalarinputs%%)
%%scalardefs%%
Cf2py intent(out) bode
      double complex zsinh, zcosh
