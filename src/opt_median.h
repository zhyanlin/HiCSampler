// hard code median search for small array
// cannot go faster unless assumptions are made

typedef double pixelvalue;
#define PIX_SORT(a, b)          \
    {                           \
        if ((a) > (b))          \
            PIX_SWAP((a), (b)); \
    }
#define PIX_SWAP(a, b)         \
    {                          \
        pixelvalue temp = (a); \
        (a) = (b);             \
        (b) = temp;            \
    }
#define PIX_MIN(a, b) ((a) > (b)) ? (b) : (a)
#define PIX_MAX(a, b) ((a) > (b)) ? (a) : (b)

/*----------------------------------------------------------------------------
    Function :   opt_med3()
    In       :   pointer to array of 3 pixel values
    Out      :   a pixelvalue
    Job      :   optimized search of the median of 3 pixel values
    Notice   :   found on sci.image.processing
    cannot go faster unless assumptions are made
    on the nature of the input signal.
    ---------------------------------------------------------------------------*/
template <class pixelvalue>
inline pixelvalue opt_med3(pixelvalue *p)
{
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[1], p[2]);
    PIX_SORT(p[0], p[1]);
    return p[1];
}

/*
    ** opt_med4()
    hacked version
    **/
template <class pixelvalue>
inline pixelvalue opt_med4(pixelvalue *p)
{
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[2], p[3]);
    PIX_SORT(p[0], p[2]);
    PIX_SORT(p[1], p[3]);
    return (p[1] + p[2]) / 2;
}

/*----------------------------------------------------------------------------
    Function :   opt_med5()
    In       :   pointer to array of 5 pixel values
    Out      :   a pixelvalue
    Job      :   optimized search of the median of 5 pixel values
    Notice   :   found on sci.image.processing
    cannot go faster unless assumptions are made
    on the nature of the input signal.
    ---------------------------------------------------------------------------*/
template <class pixelvalue>
inline pixelvalue opt_med5(pixelvalue *p)
{
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[3], p[4]);
    p[3] = PIX_MAX(p[0], p[3]);
    p[1] = PIX_MIN(p[1], p[4]);
    PIX_SORT(p[1], p[2]);
    p[2] = PIX_MIN(p[2], p[3]);
    return PIX_MAX(p[1], p[2]);
}

/*----------------------------------------------------------------------------
   Function :   opt_med6()
   In       :   pointer to array of 6 pixel values
   Out      :   a pixelvalue
   Job      :   optimized search of the median of 6 pixel values
   Notice   :   from Christoph_John@gmx.de
                based on a selection network which was proposed in
                "FAST, EFFICIENT MEDIAN FILTERS WITH EVEN LENGTH WINDOWS"
                J.P. HAVLICEK, K.A. SAKADY, G.R.KATZ
                If you need larger even length kernels check the paper
 ---------------------------------------------------------------------------*/

template <class pixelvalue>
inline pixelvalue opt_med6(pixelvalue *p)
{
    PIX_SORT(p[1], p[2]);
    PIX_SORT(p[3], p[4]);
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[2], p[3]);
    PIX_SORT(p[4], p[5]);
    PIX_SORT(p[1], p[2]);
    PIX_SORT(p[3], p[4]);
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[2], p[3]);
    PIX_SORT(p[4], p[5]);
    PIX_SORT(p[1], p[2]);
    PIX_SORT(p[3], p[4]);
    return (p[2] + p[3]) * 0.5;
    /* PIX_SORT(p[2], p[3]) results in lower median in p[2] and upper median in p[3] */
}

/*----------------------------------------------------------------------------
    Function :   opt_med7()
    In       :   pointer to array of 7 pixel values
    Out      :   a pixelvalue
    Job      :   optimized search of the median of 7 pixel values
    Notice   :   found on sci.image.processing
    cannot go faster unless assumptions are made
    on the nature of the input signal.
    ---------------------------------------------------------------------------*/
template <class pixelvalue>
inline pixelvalue opt_med7(pixelvalue *p)
{
    PIX_SORT(p[0], p[5]);
    PIX_SORT(p[0], p[3]);
    PIX_SORT(p[1], p[6]);
    PIX_SORT(p[2], p[4]);
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[3], p[5]);
    PIX_SORT(p[2], p[6]);
    p[3] = PIX_MAX(p[2], p[3]);
    p[3] = PIX_MIN(p[3], p[6]);
    p[4] = PIX_MIN(p[4], p[5]);
    PIX_SORT(p[1], p[4]);
    p[3] = PIX_MAX(p[1], p[3]);
    return PIX_MIN(p[3], p[4]);
}

template <class pixelvalue>
inline pixelvalue opt_med8(pixelvalue *p)
{
    PIX_SORT(p[0], p[4]);
    PIX_SORT(p[1], p[5]);
    PIX_SORT(p[2], p[6]);
    PIX_SORT(p[3], p[7]);
    PIX_SORT(p[0], p[2]);
    PIX_SORT(p[1], p[3]);
    PIX_SORT(p[4], p[6]);
    PIX_SORT(p[5], p[7]);
    PIX_SORT(p[2], p[4]);
    PIX_SORT(p[3], p[5]);
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[2], p[3]);
    PIX_SORT(p[4], p[5]);
    PIX_SORT(p[6], p[7]);
    PIX_SORT(p[1], p[4]);
    PIX_SORT(p[3], p[6]);
    return (p[3] + p[4]) * 0.5;
}

/*----------------------------------------------------------------------------
    Function :   opt_med9()
    In       :   pointer to an array of 9 pixelvalues
    Out      :   a pixelvalue
    Job      :   optimized search of the median of 9 pixelvalues
    Notice   :   in theory, cannot go faster without assumptions on the
    signal.
    Formula from:
    XILINX XCELL magazine, vol. 23 by John L. Smith
    The input array is modified in the process
    The result array is guaranteed to contain the median
    value
    in middle position, but other elements are NOT sorted.
    ---------------------------------------------------------------------------*/
template <class pixelvalue>
inline pixelvalue opt_med9(pixelvalue *p)
{
    PIX_SORT(p[1], p[2]);
    PIX_SORT(p[4], p[5]);
    PIX_SORT(p[7], p[8]);
    PIX_SORT(p[0], p[1]);
    PIX_SORT(p[3], p[4]);
    PIX_SORT(p[6], p[7]);
    PIX_SORT(p[1], p[2]);
    PIX_SORT(p[4], p[5]);
    PIX_SORT(p[7], p[8]);
    p[3] = PIX_MAX(p[0], p[3]);
    p[5] = PIX_MIN(p[5], p[8]);
    PIX_SORT(p[4], p[7]);
    p[6] = PIX_MAX(p[3], p[6]);
    p[4] = PIX_MAX(p[1], p[4]);
    p[2] = PIX_MIN(p[2], p[5]);
    p[4] = PIX_MIN(p[4], p[7]);
    PIX_SORT(p[4], p[2]);
    p[4] = PIX_MAX(p[6], p[4]);
    return PIX_MIN(p[4], p[2]);
}
