
#include "detectText.h"

int nbcGenerateSmart_Ex(pixel*data, int stride,int length, int max) {
    int mb_height =length;
    int mb_width =length;
    int nbc_block = 0;
    int color;
    const int nbc_diff = NBC_DIFF;
    int m, n, i, j;
    const int SIZE = 256;
    int same_value_count[257] = {0};
	pixel *p_data_tmp = data;
    int stride_diff = stride - mb_width;

    for (m = 0; m < mb_height; m++) {
        for (n = 0; n < mb_width; n++) {
            color = *p_data_tmp;
            same_value_count[color]++;
            ++p_data_tmp;
        }
        p_data_tmp += stride_diff;
    }

    same_value_count[256] = -1;
    int max_color[4] = {256,256,256,256};
    for (m = 0; m < SIZE; ++m) {
        if (same_value_count[max_color[0]] < same_value_count[m]){
            for(i = 3; i > 0; --i){
                max_color[i] = max_color[i-1];
            }
            max_color[0] = m;
        }else if (same_value_count[max_color[3]] >= same_value_count[m]){
            //do nothing
        }else{
            for(n = 0; n < 4; ++n){
                if(same_value_count[max_color[n]] < same_value_count[m]){
                    for(i = 3; i > n; --i){
                        max_color[i] = max_color[i-1];
                    }
                    max_color[n] = m;
                    break;
                }
            }
        }
    }
    nbc_block = 0;
    for(i = 0; i < 4; ++i){
        nbc_block += same_value_count[max_color[i]];
        same_value_count[max_color[i]] = 0;
        for(j = 1; j < nbc_diff; ++j){
            if(max_color[i]-j >= 0){
                nbc_block += same_value_count[max_color[i]-j];
                same_value_count[max_color[i]-j] = 0;
            }
            if(max_color[i]+j < 256){
                nbc_block += same_value_count[max_color[i]+j];
                same_value_count[max_color[i]+j] = 0;
            }
        }
    }
    return nbc_block;
}

int nhgGenerateSmart_Ex(pixel *data, int stride,int length,int max_value) {
    int mb_height = length;
    int mb_width = length;
    int i, j;
    int nhg_block = 0;
    int grads = 0;
	pixel *p_data_tmp = data;
    int stride_diff = stride - mb_width;
    int mb_width_last_p = mb_width - 1;

    for (i = 1; i < mb_height-1; ++i) {
        for (j = 1; j < mb_width-1; ++j,++p_data_tmp) {
            grads = abs(*p_data_tmp - *(p_data_tmp +1));
            if (grads > NGH_DIFF) {
                ++nhg_block;
                continue;
            }
            grads = abs(*p_data_tmp - *(p_data_tmp + stride));
            if (grads > NGH_DIFF) {
                ++nhg_block;
                continue;
            }
            grads = abs(*p_data_tmp - *(p_data_tmp - stride));
            if (grads > NGH_DIFF) {
                ++nhg_block;
                continue;
            }
            grads = abs(*p_data_tmp - *(p_data_tmp - 1));
            if (grads > NGH_DIFF) {
                ++nhg_block;
                continue;
            }
        }
        p_data_tmp += stride_diff;
    }
    return nhg_block;
}
int decideTextBlock(pixel* p_src, int stride_p_src, int length){
    int mb_nhg_count = 0;
    int mb_nbc_count = 0;
    int ret_val = 0;
    mb_nbc_count = nbcGenerateSmart_Ex(p_src, stride_p_src,length,NBC_TH0);
    mb_nhg_count = nhgGenerateSmart_Ex(p_src, stride_p_src,length,NGH_TH);
    if (mb_nbc_count > 160) {
        ret_val = 1;
        if (mb_nhg_count < 40) {
            ret_val = 0;
        }
    }else {
        ret_val = 1;
        if (mb_nhg_count < 50) {
            ret_val = 0;
        }
    }
    return ret_val;
}


