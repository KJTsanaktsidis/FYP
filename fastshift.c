int fast_pad_shift(double* y1, double* y2, int len)
{
    //shift y2 on y1

    double min_lsq_sum = -2;
    int best_shift = 0;

    //want to shift by values of dx from -len to len
    for (int dx = -len; dx <= len; dx++)
    {
        double lsq_sum = 0;

        //each pair
        for (int i = 0; i < len; i++)
        {
            double y1val = y1[i];
            int y2index = i + dx;
            double y2val;

            if (y2index < 0)
                y2val = 1;
            else if (y2index >= len)
                y2val = 0;
            else
                y2val = y2[y2index];

            lsq_sum += (y2val - y1val) * (y2val - y1val);
        }

        if (lsq_sum < min_lsq_sum || min_lsq_sum < -1)
        {
            min_lsq_sum = lsq_sum;
            best_shift = dx;
        }
    }

    //Make the convention that left-shifts are negative, right-shifts are positive
    return -best_shift;
}
