/*.BA*/



/*.FE{C 4.15.1}{Fehler und Kondition}{Fehler und Kondition}*/

/*.BE*/
/* -------------------------- MODUL fcond.c ------------------------- */

#include "basis.h"
#include "vmblock.h"
#include "u_proto.h"

/*.BA*/

/*.BE*/
REAL hcond              /* Hadamardsche Konditionszahl ...............*/
/*.BA*/
/*.IX{hcond}*/
/*.BE*/
             (
              int     n,          /* Dimension der Matrix ............*/
              REAL *  mat[]       /* Eingabematrix ...................*/
             )
/*.BA*/

/*====================================================================*
 *                                                                    *
 *  hcond bestimmt die Hadamardsche Konditionszahl einer n x n        *
 *  Matrix.                                                           *
.BE*)
 *  Ist der Rueckgabewert von hcond() sehr viel kleiner als 1, so ist *
 *  die Matrix schlecht konditioniert. Die Loesung eines linearen     *
 *  Gleichungssystems wird dann ungenau.                              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension von mat.                                   *
 *      mat      REAL   *mat[n];                                      *
 *               n x n Matrix, deren Konditionszahl zu bestimmen ist. *
 *                                                                    *
 *   Rueckgabewert:                                                   *
 *   =============                                                    *
 *      REAL     < 0.0: Fehler                                        *
 *                 = -1.0 :  n < 1                                    *
 *                 = -2.0 :  zu wenig Speicher                        *
 *                 = -3.0 :  Matrix ist singulaer (det = 0.0)         *
 *                                                                    *
 *               >= 0.0 :                                             *
 *               Hadamardsche Konditionszahl von mat.                 *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Funktionen:                                             *
 *   ===================                                              *
 *                                                                    *
 *      int gaudec ():    Zerlegung von mat in LU-Form.               *
 *      void *vmalloc (): Speicher fuer Vektor bzw. Matrix anfordern. *
 *      void vmfree ():   Liste von Vektoren und Matrizen freigeben.  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Konstanten: NULL                                        *
 *   ===================                                              *
 *                                                                    *
 *   Macros: ABS, SQRT                                                *
 *   ======                                                           *
.BA*)
 *====================================================================*/
/*.BE*/
{
  register int j, i;
  REAL     temp, cond, **lu;
  int      rc, signd, *perm;
  void *vmblock;

  if (n < 1) return (-ONE);

                                          /* Speicher fuer die        */
  vmblock = vminit();                     /* Gausszerlegung anfordern */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (-TWO);
  }

  rc = gaudec (n, mat, lu, perm, &signd);       /* Zerlegung in lu    */
  if (rc != 0 || signd == 0)                    /* berechnen          */
  {
    vmfree(vmblock);
    return (-THREE);
  }

  cond = ONE;                                   /* Konditionszahl     */
  for (i = 0; i < n; i++)                       /* bestimmen          */
  {
    for (temp = ZERO, j = 0; j < n; j++)
      temp += SQR (mat[i][j]);
    cond *= lu[i][i] / SQRT (temp);
  }

  vmfree(vmblock);                              /* Speicher freigeben */

  return (ABS (cond));
}
/*.BA*/



/*.FE{C 4.15.2}{Konditionssch"atzung}{Konditionssch"atzung}*/

/*.FE{}{Konditionssch"atzung nach Cline}
       {Konditionssch"atzung nach Cline}*/

/*.BE*/
REAL ccond              /* Konditionsschaetzung nach Cline ...........*/
/*.BA*/
/*.IX{ccond}*/
/*.BE*/
             (
              int     n,          /* Dimension der Matrix ............*/
              REAL *  mat[]       /* Eingabematrix ...................*/
             )
/*.BA*/

/*====================================================================*
 *                                                                    *
 *  ccond bestimmt eine Schaetzung der Konditionszahl cond (mat) einer*
 *  n x n Matrix mat nach dem Verfahren von Cline.                    *
.BE*)
 *                               -1                                   *
 *  Es ist cond (A) = | A | * | A  |, wobei hier | | fuer die         *
 *  Maximumnorm steht.                                                *
 *                                                                    *
 *  Ein grosser Wert fuer cond (A) deutet schlechte Kondition der     *
 *  Matrix A an. Loesungen von Gleichungssystemen werden i.a. ungenau.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension von mat.                                   *
 *      mat      REAL   *mat[n];                                      *
 *               n x n Matrix, deren Konditionszahl zu schaetzen ist. *
 *                                                                    *
 *   Rueckgabewert:                                                   *
 *   =============                                                    *
 *      REAL     < 0.0: Fehler                                        *
 *                 = -1.0 :  n < 1                                    *
 *                 = -2.0 :  zu wenig Speicher                        *
 *                 = -3.0 :  Matrix ist singulaer (det = 0.0)         *
 *                                                                    *
 *               >= 0.0 :                                             *
 *               Cline's Schaetzung der Konditionszahl von mat.       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Funktionen:                                             *
 *   ===================                                              *
 *                                                                    *
 *      int gaudec ():    Zerlegung von mat in LU-Form.               *
 *      void *vmalloc (): Speicher fuer Vektor bzw. Matrix anfordern. *
 *      void vmfree ():   Liste von Vektoren und Matrizen freigeben.  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Konstanten: NULL                                        *
 *   ===================                                              *
 *                                                                    *
 *   Macros: ABS                                                      *
 *   ======                                                           *
.BA*)
 *====================================================================*/
/*.BE*/
{
  register int i, j, k;
  REAL     **lu, *x, *y, *z,
           v, smi, spl, sum, xnorm, znorm, matnorm;
  int      rc, signd, *perm;
  void *vmblock;

  if (n < 1) return (-ONE);
                                         /* Speicher fuer die         */
  vmblock = vminit();                    /* Gausszerlegung allokieren */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (-TWO);
  }

  rc = gaudec (n, mat, lu, perm, &signd);       /* Zerlegung in lu    */
  if (rc != 0 || signd == 0)                    /* berechnen          */
  {
    vmfree(vmblock);                            /* Speicher freigeben */
    return (-THREE);
  }

  x = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  y = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  z = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))                    /* kein Speicher .....*/
  {
    vmfree(vmblock);
    return (-THREE);
  }
                                    /* Bestimme x = (+-1,+-1...,+-1)  */
  x[0] = ONE;                       /* so, dass y "maximal" wird      */
  y[0] = ONE / lu[0][0];
  for (i = 1; i < n; i++)
    y[i] = - lu[0][i] * y[0] / lu[i][i];

  for (k = 1; k < n; k++)
  {
    v = ONE / lu[k][k];
    x[k] = y[k] - v;
    y[k] += v;
    smi = ABS (x[k]);
    spl = ABS (y[k]);
    for (i = k + 1; i < n; i++)
    {
      v = lu[k][i] / lu[i][i];
      x[i] = y[i] - v * x[k];
      y[i] -= v * y[k];
      smi += ABS (x[i]);
      spl += ABS (y[i]);
    }

    if (smi > spl)
    {
      for (i = k; i < n; i++) y[i] = x[i];
      x[k] = -ONE;
    }
    else
      x[k] = ONE;
  }

  for (i = n - 1; i >= 0; i--)         /* Rueckwaertselimination .....*/
  {
    z[i] = y[i];
    for (j = i + 1; j < n; j++)
      z[i] -= lu[j][i] * y[j];
  }

  znorm = ZERO;                        /* Normen bestimmen ...........*/
  xnorm = ZERO;
  matnorm = ZERO;
  for (i = 0; i < n; i++)
  {
    if (ABS (z[i]) > znorm)            /* Maximumnorm von z ..........*/
      znorm = ABS (z[i]);
    if (ABS (x[i]) > xnorm)            /* Maximumnorm von x ..........*/
      xnorm = ABS (x[i]);

    sum = ZERO;
    for (j = 0; j < n; j++)            /* Maximumnorm von mat ........*/
      sum += ABS (mat[i][j]);
    if (sum > matnorm)
      matnorm = sum;
  }

  vmfree(vmblock);                              /* Speicher freigeben */

  return (matnorm * (znorm / xnorm));
}


/*.BA*/



/*.FE{}{Konditionssch"atzung nach Forsythe/Moler}
       {Konditionssch"atzung nach Forsythe/Moler}*/

/*.BE*/
REAL fcond              /* Konditionsschaetzung nach Forsythe/Moler ..*/
/*.BA*/
/*.IX{fcond}*/
/*.BE*/
             (
              int     n,          /* Dimension der Matrix ............*/
              REAL *  mat[]       /* Eingabematrix ...................*/
             )
/*.BA*/

/*====================================================================*
 *                                                                    *
 *  fcond bestimmt eine Schaetzung der Konditionszahl cond (mat) einer*
 *  n x n Matrix mat nach dem Verfahren von Forsythe/Moler.           *
.BE*)
 *                               -1                                   *
 *  Es ist cond (A) = | A | * | A  |, wobei hier | | fuer die         *
 *  Maximumnorm steht.                                                *
 *                                                                    *
 *  Ein grosser Wert fuer cond (A) deutet schlechte Kondition der     *
 *  Matrix A an. Loesungen von Gleichungssystemen werden i.a. ungenau.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension von mat.                                   *
 *      mat      REAL   *mat[n];                                      *
 *               n x n Matrix, deren Konditionszahl zu schaetzen ist. *
 *                                                                    *
 *   Rueckgabewert:                                                   *
 *   =============                                                    *
 *      REAL     < 0.0: Fehler                                        *
 *                 = -1.0 :  n < 1                                    *
 *                 = -2.0 :  zu wenig Speicher                        *
 *                 = -3.0 :  Matrix ist singulaer (det = 0.0)         *
 *                                                                    *
 *               >= 0.0 :                                             *
 *               Schaetzung der Konditionszahl von mat nach           *
 *               Forsythe/Moler.                                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Funktionen:                                             *
 *   ===================                                              *
 *                                                                    *
 *      int gauss () :    Loesung eines Gleichungssystems.            *
 *      void *vmalloc (): Speicher fuer Vektor bzw. Matrix anfordern. *
 *      void vmfree ():   Liste von Vektoren und Matrizen freigeben.  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Konstanten: NULL                                        *
 *   ===================                                              *
 *                                                                    *
 *   Macros: ABS, MACH_EPS                                            *
 *   ======                                                           *
.BA*)
 *====================================================================*/
/*.BE*/
{
  register int i, j;
  REAL     **lu, *x, *b, *r, nom, denom;
  int      rc, signd, *perm;
  LONG_REAL sum;
  void *vmblock;

  if (n < 1) return (-ONE);
                                          /* Speicher fuer die        */
  vmblock = vminit();                     /* Gausszerlegung anfordern */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));
  x = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  b = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  r = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (-TWO);
  }

  for (i = 0; i < n; i++) b[i] = ONE;      /* b = 1-Vektor            */

  /* Loese mat * x = b ...............................................*/
  rc = gauss (0, n, mat, lu, perm, b, x, &signd);

  if (rc)                                  /* Matrix singulaer        */
  {
    vmfree(vmblock);
    return (-THREE);
  }

  for (i = 0; i < n; i++)                  /* Residuen bestimmen      */
  {                                        /* mit long double         */
    sum = (LONG_REAL) b[i];
    for (j = 0; j < n; j++)
      sum -= (LONG_REAL) mat[i][j] * (LONG_REAL) x[j];
    r[i] = (REAL) sum;
  }

  /* Loese mat * b = r, d.h. eine Nachiteration ......................*/
  rc = gauss (2, n, mat, lu, perm, r, b, &signd);

  if (rc != 0 || signd == 0)               /* Sollte nie passieren    */
  {
    vmfree(vmblock);
    return (-THREE);
  }

  denom = nom = ZERO;                      /* Max-Normen bestimmen    */
  for (i = 0; i < n; i++)
  {
    if (denom < ABS (b[i])) denom = ABS (b[i]);
    if (nom   < ABS (x[i])) nom   = ABS (x[i]);
  }

  vmfree(vmblock);

  return (denom / nom / MACH_EPS);
}

#define MAXITER 30      /* Maximalzahl von Nachiterationen ...........*/
/*.BA*/



/*.FE{C 4.15.4}{Nachiteration}{Nachiteration}*/

/*.BE*/
int gausoli              /* Gauss Loesung ............................*/
/*.BA*/
/*.IX{gausoli}*/
/*.BE*/
            (
             int     n,            /* Dimension der Matrix ...........*/
             REAL *  mat[],        /* Ausgangsmatrix .................*/
             REAL *  lumat[],      /* Eingabematrix (LU) .............*/
             int     perm[],       /* Zeilenvertauschungen ...........*/
             REAL    b[],          /* Rechte Seite ...................*/
             REAL    x[]           /* Loesung ........................*/
            )
/*.BA*/

/*====================================================================*
 *                                                                    *
 *  gausoli bestimmt die Loesung x des linearen Gleichungssystems     *
 *  lumat * x = b mit der n x n Koeffizientenmatrix lumat, wobei      *
 *  lumat in zerlegter Form ( LU - Dekomposition ) vorliegt, wie      *
 *  sie von gaudec als Ausgabe geliefert wird.                        *
 *  Die Loesung wird durch Nachiteration verbessert.                  *
.BE*)
 *  Die Nachiteration wird abgebrochen, wenn die relative Ver-        *
 *  besserung kleiner als 2*MACH_EPS ist oder die Norm der Re-        *
 *  siduen ansteigt oder die maximale Iterationszahl erreicht ist.    *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension von lumat,                                 *
 *               Anzahl der Komponenten des b-Vektors, des Loe-       *
 *               sungsvektors x, des Permutationsvektors perm.        *
 *      mat      REAL   *mat[];                                       *
 *               Matrix des Gleichungssystems. Diese wird als Vektor  *
 *               von Zeigern uebergeben.                              *
 *      lumat    REAL   *lumat[];                                     *
 *               LU-Dekompositionsmatrix, wie sie von gaudec          *
 *               geliefert wird.                                      *
 *               Achtung: mat und lumat muessen hierbei verschieden   *
 *                        gewahlt sein !!!                            *
 *      perm     int perm[];                                          *
 *               Permutationsvektor, der die Zeilenvertauschungen     *
 *               von lumat enthaelt.                                  *
 *      b        REAL   b[];                                          *
 *               Rechte Seite des Gleichungssystems.                  *
 *                                                                    *
 *   Ausgabeparameter:                                                *
 *   ================                                                 *
 *      x        REAL   x[];                                          *
 *               Loesungsvektor des Systems.                          *
 *                                                                    *
 *   Rueckgabewert:                                                   *
 *   =============                                                    *
 *      =-1      Maximale Nachiterationszahl (MAXITER) erreicht       *
 *      = 0      alles ok                                             *
 *      = 1      n < 1 gewaehlt oder unzulaessige Eingabeparameter    *
 *      = 2      zu wenig Speicherplatz                               *
 *      = 3      unzulaessige Zerlegungsmatrix                        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Funktionen:                                             *
 *   ===================                                              *
 *                                                                    *
 *      int gausol ():    Ausgangsloesung bestimmen                   *
 *      void *vmalloc (): Speicher fuer Vektor bzw. Matrix anfordern. *
 *      void vmfree ():   Liste von Vektoren und Matrizen freigeben.  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Benutzte Konstanten: NULL, MAXROOT, MACH_EPS, MAXITER            *
 *   ===================                                              *
 *                                                                    *
.BA*)
 *====================================================================*/
/*.BE*/
{
  int       i, j, k, rc;
  REAL      *r, *z, maxx, maxz, oldmaxz, eps;
  LONG_REAL sumld;
  void *vmblock;

  if (n < 1) return (1);
  if (mat == lumat) return (1);
                                 /* Loesen des Systems mit gauss .....*/
  if ((rc = gausol (n, lumat, perm, b, x)) != 0)
    return rc;

  eps = (REAL) (TWO * MACH_EPS);
  oldmaxz = MAXROOT;

  vmblock = vminit();
  z = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  r = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))
    return 2;

  for (k = 1; ; k++)
  {
    maxx = ZERO;
    for (i = 0; i < n; i++)      /* Residuen "genauer" bestimmen .....*/
    {
      sumld = (LONG_REAL) b[i];
      for (j = 0; j < n; j++)
        sumld -= (LONG_REAL) mat[i][j] * (LONG_REAL) x[j];

      r[i] = (REAL) sumld;
      if (ABS (x[i]) > maxx) maxx = ABS (x[i]);
    }

    rc = gausol (n, lumat, perm, r, z);    /* Loese mat * z = r ......*/
    if (rc) break;

    maxz = ZERO;          /* x korrigieren, max (ABS(z[i])) bestimmen */
    for (i = 0; i < n; i++)
    {
      x[i] += z[i];
      if (ABS (z[i]) > maxz) maxz = ABS (z[i]);
    }

    if (maxz < eps * maxx)    /* Ende pruefen ........................*/
    {
      rc = 0;
      break;
    }

    if (oldmaxz < maxz)       /* Divergenz ? .........................*/
    {
      rc = 0;
      break;
    }

    if (k >= MAXITER)         /* Maximale Iterationszahl erreicht ....*/
    {
      rc = -1;
      break;
    }

    oldmaxz = maxz;           /* Letzte Maximumnorm von z merken .....*/
  }  /* end of k */

  vmfree(vmblock);            /* Allokierten Speicher freigeben ......*/

  return rc;
}

/* -------------------------- ENDE fcond.c -------------------------- */
