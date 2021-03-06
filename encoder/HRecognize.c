//
// Created by xiangyy on 2017/9/2.
//
#include "encodeText.h"
#include "HRecognize.h"


//限定Tile的宽度和高度必须是MB_WIDTH和MB_HEIGHT的整数倍
//按照以上的限制重置Tile的尺寸
void resetTileSize(int frame_width, int frame_height, int *tile_width, int *tile_height){
    int aligned_frame_width = (((frame_width-1)>>BIT_SHIFT_MB_WIDTH)+1)*MB_WIDTH;
    int aligned_frame_height = (((frame_height-1)>>BIT_SHIFT_MB_HEIGHT)+1)*MB_HEIGHT;
    *tile_width = (((*tile_width-1)>>BIT_SHIFT_MB_WIDTH)+1)*MB_WIDTH;
    *tile_height = (((*tile_height-1)>>BIT_SHIFT_MB_HEIGHT)+1)*MB_HEIGHT;
    if(*tile_width>aligned_frame_width){
        *tile_width = aligned_frame_width;
    }
    if(*tile_height>aligned_frame_height){
        *tile_height = aligned_frame_height;
    }
}

//分割Tile的方法
//根据Frame和Tile的尺寸，计算共有多少个Tile，计算每个Tile的尺寸和位置
int generateTileSizeInfoTable(int frame_width, int frame_height, int tile_width, int tile_height, struct TileSizeAndPosition **table) {
    int tile_count = 0; //tile的个数
    int ii=0, jj=0;
    int tile_in_row = (frame_width-1)/tile_width + 1;
    int tile_in_column = (frame_height-1)/tile_height + 1;
    tile_count = tile_in_row * tile_in_column;
    *table = (struct TileSizeAndPosition *)malloc(tile_count*sizeof(struct TileSizeAndPosition));
    memset(*table, 0, tile_count*sizeof(struct TileSizeAndPosition));
    for(jj=0; jj<tile_in_column;jj++) {
        for(ii=0; ii<tile_in_row;ii++) {
            struct TileSizeAndPosition *pt = *table + jj*tile_in_row + ii;
            pt->tile_id = jj*tile_in_row + ii;
            pt->pos_x = ii*tile_width;
            pt->pos_y = jj*tile_height;
            if(ii != (tile_in_row - 1)){
                //不是最右边的Tile
                pt->width = tile_width;
            }else{
                pt->width = frame_width - ii * tile_width;
            }
            if (jj != (tile_in_column - 1)) {
                //不是最右边的Tile
                pt->height = tile_height;
            } else {
                pt->height = frame_height - jj * tile_height;
            }
            pt->stride = frame_width + 2*IMAGE_PADDING_SIZE;
            pt->tile_in_row = tile_in_row;
            pt->tile_in_column = tile_in_column;
            pt->frame_width = frame_width;
            pt->frame_height = frame_height;
            pt->tile_width = tile_width;
            pt->tile_height = tile_height;
            pt->mb_count_in_tile_width = (pt->width-1) / MB_WIDTH + 1;
            pt->mb_count_in_tile_height = (pt->height-1) / MB_HEIGHT + 1;
            pt->tile_index_x = ii;
            pt->tile_index_y = jj;
        }
    }
    return tile_count;
}

void freeTileSizeInfoTable(struct TileSizeAndPosition *table) {
    safe_free(table);
    return ;
}

int generateTilesFromFrame(struct RecognizeContext *const context) {
    int success = -1;
    struct TileSizeAndPosition *tile_info_table = context->tile_info_table;
    int tile_count = context->tile_count;
    if((tile_count<1)||(NULL == tile_info_table)){
        //计算每个Tile的尺寸没有成功，退出
        return success;
    }
    struct Tile *tile_array = (struct Tile *)calloc(sizeof(struct Tile), (size_t)tile_count);
    context->tile_array_of_current_frame = tile_array;
    for(int i=0;i<tile_count;i++){
        //逐个Tile完成Tile的尺寸的赋值
        tile_array[i].frame_id = -1; //置为缺省值
        tile_array[i].tile_id = i;
        tile_array[i].width = tile_info_table[i].width;  //有效像素宽度
        tile_array[i].height = tile_info_table[i].height; //有效像素高度
        tile_array[i].pos_x = tile_info_table[i].pos_x; //在Frame中的x位置。
        tile_array[i].pos_y = tile_info_table[i].pos_y; //在Frame中的y位置。
        tile_array[i].stride = tile_info_table[i].stride;
        tile_array[i].offset = IMAGE_PADDING_SIZE * tile_array[i].stride + IMAGE_PADDING_SIZE; //图像指针相对于延拓后图像的偏移量
    }
    return success;
}


//为RecognizeContext结构体的内部指针分配内存空间
int allocMemoryForRecognizeContext(const int frame_width, const int frame_height,
                                   const int tile_width, const int tile_height,
                                   struct RecognizeContext* const context){
    int success = -1;
    //识别视频的策略
    //初始化context的mb_mode
    context->frame_width = frame_width;
    context->frame_height = frame_height;
    context->tile_width = tile_width;
    context->tile_height = tile_height;
    //重置Tile的尺寸保证Tile的宽和高是MB_WIDTH和MB_HEIGHT的整数倍
    resetTileSize(frame_width, frame_height, &(context->tile_width), &(context->tile_height));

    //得到Tile的数量
    context->tile_count = generateTileSizeInfoTable(context->frame_width, context->frame_height, \
                                                    context->tile_width, context->tile_height, \
                                                    &(context->tile_info_table));

    //产生Tile的数组
    success = generateTilesFromFrame(context);

    int mb_in_row = ((frame_width-1)/MB_WIDTH + 1);
    int mb_in_column = ((frame_height-1)/MB_HEIGHT + 1);
    int mb_count = mb_in_row * mb_in_column;
    int mb_stride = mb_in_row + 2*MB_MODE_PADDING_WIDTH;
    context->mb_in_row = mb_in_row;
    context->mb_in_column = mb_in_column;
    context->mb_count = mb_in_row * mb_in_column;
    context->mb_stride = mb_stride;
    context->mb_offset = MB_MODE_PADDING_WIDTH*mb_stride + MB_MODE_PADDING_WIDTH;
    int mb_count_padding = mb_stride*(mb_in_column+2*MB_MODE_PADDING_WIDTH);

    //分配特征索引表
    context->feature_meta_in_one_frame = (FeatureMeta *)calloc(sizeof(FeatureMeta), 1);
    context->feature_meta_in_another_frame = (FeatureMeta *)calloc(sizeof(FeatureMeta), 1);
    //要将特征数据的帧号修改为-1，否则，第0帧的帧号就与初始值相同，导致错误
    context->frame_id_of_one_frame = -1;
    context->frame_id_of_another_frame = -1;
    context->mbtype =(enum MBType *)calloc(sizeof(enum MBType), (size_t)mb_count_padding);
    context->mbtype += context->mb_offset;
    context->feature_meta = (FeatureMeta *)calloc(sizeof(FeatureMeta), 1);
    return success;
}

//释放为RecognizeContext结构体的内部指针分配内存空间
void freeMemoryForRecognizeContext(struct RecognizeContext* const context){
    //释放Tile的划分信息
    safe_free(context->tile_info_table);
    safe_free(context->tile_array_of_current_frame);
    //释放特征索引表
    safe_free(context->feature_meta_in_one_frame);
    safe_free(context->feature_meta_in_another_frame);
    context->mbtype -= context->mb_offset;
    safe_free(context->mbtype);
    safe_free(context->feature_meta);
    return;
}

//初始化识别模块的上下文
int initialRecognizeContext(const int frame_width, const int frame_height, \
                            const int tile_width, const int tile_height, \
                            enum RecognizeNoiseThreshold noise_threshold, \
                            struct RecognizeContext **context){
    int success = 0;
    //将后续使用的内存一次性分配完毕，以后可以重复使用
    //开始检测一下Tile和Frame的尺寸是否产生了变化，如果变化了，释放已经分配的内存，然后重新分配
    struct RecognizeContext *recognize_context = (struct RecognizeContext *)calloc(sizeof(struct RecognizeContext), 1);
    if(NULL==recognize_context){
        printf("%s:%d malloc recognize_context failed.",__FILE__, __LINE__);
        return -1;
    }
    *context = recognize_context;
    success = allocMemoryForRecognizeContext(frame_width, frame_height, tile_width, tile_height, recognize_context);
    recognize_context->noise_threshold = noise_threshold;

	for (int i = 0; i<recognize_context->mb_in_column; i++) {
		memset(recognize_context->mbtype + i*recognize_context->mb_stride, DEFAULT, (recognize_context->mb_in_row) * sizeof(enum MBType));
	}
    return 0;
}

//从YUV文件中读取特定的一帧内容到frame中，其中frame的大小已经指定好，内存已经提前申请，内存大小是按照存放Padding后的图像
int readFrameFromYUVFile(const char *filename, struct Frame *frame, int frame_no) {
    int success = -1;
    unsigned char *image_y = frame->data_component1;
    unsigned char *image_u = frame->data_component2;
    unsigned char *image_v = frame->data_component3;
    int frame_width = frame->width;
    int frame_height = frame->height;
    int frame_stride = frame->stride;
    FILE * pFile;
    long skip_size = frame_no*frame_width*frame_height*3;
    long read_size = frame_width;
    size_t result;
    int i;

    pFile = fopen (filename , "rb");
    if (pFile == NULL)
    {
        printf("Error opening file |%s|. \n", filename);
        success = -1;
    }
    else
    {
        if (IS_YUV420P) {
            int frame_height2 = (frame_height + 1) / 2;
            int frame_width2 = (frame_width + 1) / 2;
            int frame_stride2 = (frame_stride + 1) / 2;
            long read_size2 = frame_width2;

            long skip_size2 = frame_no * (frame_width * frame_height + frame_height2 * frame_width2 * 2);
            fseek(pFile, skip_size2, SEEK_SET);
            //y component
            for (i = 0; i < frame_height; i++) {
                //每次读一行，读完一行加stride
                result = fread(image_y, 1, (size_t) read_size, pFile);
                image_y += frame_stride;
                if ((result != read_size)) {
                    printf("Reading error.\n");
                    success = -1;
                } else {
                    success = 0;
                }
            }
            // u component
            for (i = 0; i < frame_height2; i++) {
                //每次读一行，读完一行加stride
                result = fread(image_u, 1, (size_t) read_size2, pFile);
                image_u += frame_stride2;
                if ((result != read_size2)) {
                    printf("Reading error.\n");
                    success = -1;
                } else {
                    success = 0;
                }
            }
            // v component
            for (i = 0; i < frame_height2; i++) {
                //每次读一行，读完一行加stride
                result = fread(image_v, 1, (size_t) read_size2, pFile);
                image_v += frame_stride2;
                if ((result != read_size2)) {
                    printf("Reading error.\n");
                    success = -1;
                } else {
                    success = 0;
                }
            }
        } else {
            fseek(pFile , skip_size , SEEK_SET);
            for (i = 0; i < frame_height; i++) {
                //每次读一行，读完一行加stride
                result = fread (image_y, 1, (size_t)read_size, pFile);
                image_y += frame_stride;
                if ((result != read_size)) {
                    printf("Reading error.\n");
                    success = -1;
                } else {
                    success = 0;
                }
            }

            for (i = 0; i < frame_height; i++) {
                //每次读一行，读完一行加stride
                result = fread (image_u, 1, (size_t)read_size, pFile);
                image_u += frame_stride;
                if ((result != read_size)) {
                    printf("Reading error.\n");
                    success = -1;
                } else {
                    success = 0;
                }
            }
            for (i = 0; i < frame_height; i++) {
                //每次读一行，读完一行加stride
                result = fread (image_v, 1, (size_t)read_size, pFile);
                image_v += frame_stride;
                if ((result != read_size)) {
                    printf("Reading error.\n");
                    success = -1;
                } else {
                    success = 0;
                }
            }
            //printf("Reading one tile data from YUV file success.\n");

        }
    }
    if(NULL!=pFile){
        fclose (pFile);
    }
    return success;
}

//设置Frame的数据结构并分配内存
void allocAndSetFrame(int frame_width, int frame_height, int frame_id, struct Frame **pp_frame) {
    struct Frame *p_frame = (struct Frame *)calloc(sizeof(struct Frame), 1);
    if(NULL==p_frame){
        printf("%s:%d allocAndSetFrame(): malloc failed.\n", __FILE__, __LINE__);
        return;
    }
    *pp_frame = p_frame;
    p_frame->frame_id = frame_id;
    p_frame->width = frame_width;
    p_frame->height = frame_height;
    p_frame->stride = frame_width + 2 * IMAGE_PADDING_SIZE;
    p_frame->component_memory_size = (frame_width+2*IMAGE_PADDING_SIZE)*(frame_height+2*IMAGE_PADDING_SIZE);
    p_frame->offset = IMAGE_PADDING_SIZE*p_frame->stride + IMAGE_PADDING_SIZE;

    //分配空间
    p_frame->data_component1 = (unsigned char *)calloc(1, (size_t)p_frame->component_memory_size);

    if(GLOBAL_COLOR_FORMAT == kYUV444){
        p_frame->frame_format = kYUV444;
        p_frame->stride2 = p_frame->stride;
        p_frame->offset2 = p_frame->offset;
        p_frame->data_component2 = (unsigned char *)calloc(1, (size_t)p_frame->component_memory_size);
        p_frame->data_component3 = (unsigned char *)calloc(1, (size_t)p_frame->component_memory_size);
    }else{
        p_frame->frame_format = kYUV420p;
        p_frame->stride2 = (frame_width+1)/2 + 2 * IMAGE_PADDING_SIZE2;
        p_frame->offset2 = IMAGE_PADDING_SIZE2*p_frame->stride2 + IMAGE_PADDING_SIZE2;
        int component_memory_size2 = ((frame_width+1)/2+2*IMAGE_PADDING_SIZE2)*((frame_height+1)/2+2*IMAGE_PADDING_SIZE2);
        p_frame->data_component2 = (unsigned char *)calloc(1, (size_t)component_memory_size2);
        p_frame->data_component3 = (unsigned char *)calloc(1, (size_t)component_memory_size2);
    }

    if(NULL==p_frame->data_component1||
       NULL==p_frame->data_component2||
       NULL==p_frame->data_component3){
        printf("%s:%d allocAndSetFrame(): malloc failed.\n", __FILE__, __LINE__);

        safe_free(p_frame->data_component1);
        safe_free(p_frame->data_component2);
        safe_free(p_frame->data_component3);

        safe_free(p_frame);
        return;
    }
    p_frame->data_component1 += p_frame->offset;
    p_frame->data_component2 += p_frame->offset2;
    p_frame->data_component3 += p_frame->offset2;
    return;
}


void freeRecognizeContext(struct RecognizeContext *context){
    freeMemoryForRecognizeContext(context);
    safe_free(context);
    return;
}

void freeFrame(struct Frame *frame, bool need_free) {
    if (frame == NULL) {
        return;
    }
    if (need_free) {
        if (frame->offset) {
            unsigned char *ptr = frame->data_component1 - frame->offset;
            safe_free(ptr);
            ptr = frame->data_component2 - frame->offset2;
            safe_free(ptr);
            ptr = frame->data_component3 - frame->offset2;
            safe_free(ptr);
        } else {
            safe_free(frame->data_component1);
            safe_free(frame->data_component2);
            safe_free(frame->data_component3);
        }
    }
    safe_free(frame);
}

void recordFrameIntoYUVFile(char *filename, struct Frame *frame){
    //记录图像到文件中
    FILE * pRecordFile;
    int ii;
    pRecordFile = fopen (filename, "ab+");
    if (pRecordFile == NULL)
    {
        printf("Error opening RECORD file |%s|\n.", filename);
    }
    else
    {
        fseek(pRecordFile, 0, SEEK_END);
        unsigned char *p = frame->data_component1;
        for (ii=0;ii<frame->height;ii++)
        {
            fwrite(p, 1, frame->width, pRecordFile);
            p += frame->stride;
        }
        if(IS_YUV420P){
            p = frame->data_component2;
            int height2 = (frame->height+1)/2;
            int width2  = (frame->width+1)/2;
            int stride2 = (frame->stride+1)/2;
            for (ii=0;ii<height2;ii++)
            {
                fwrite(p, 1, width2, pRecordFile);
                p += stride2;
            }

            p = frame->data_component3;
            for (ii=0;ii<height2;ii++)
            {
                fwrite(p, 1, width2, pRecordFile);
                p += stride2;
            }
        }else{
            p = frame->data_component2;
            for (ii=0;ii<frame->height;ii++)
            {
                fwrite(p, 1, frame->width, pRecordFile);
                p += frame->stride;
            }

            p = frame->data_component3;
            for (ii=0;ii<frame->height;ii++)
            {
                fwrite(p, 1, frame->width, pRecordFile);
                p += frame->stride;
            }
        }

        fclose(pRecordFile);
    }
    //printf("Record a frame over!\n");
}

static unsigned short calcHashKeyByStepXor13(unsigned char *message, unsigned int len)
{
    unsigned int step = 5;
    unsigned short sum = 0;

    for (int i = 0; i < len; i++) {
        sum += message[i];
    }
    unsigned short average = sum / len;
    unsigned short value_key = average;
    unsigned short value_mask = 0;

    for (int i = 0; i < step; i++) {
        unsigned short value = 0;
        for (int j = i; j < len; j = j + step){
            value = value ^ (message[j] >= average);
        }
        value_mask = (value_mask<<1) | (value & 0x1);
    }
    value_key = (value_key << 8) | value_mask;

    return value_key;
}

int findAllFeaturesInFrame(struct RecognizeContext * context, const struct Frame *frame, int *feature_count,
                           int *collision_count, FeatureTable *feature_table, FeatureIndexTable *index_table){
    int y,x,this_pos,top_line;
    unsigned char  * data = frame->data_component1;
    const int stride = frame->stride;
    FeatureItem *feature_start = context->feature_meta->feature_table.feature_item;
    FeatureItem * feature_item;
    unsigned char  *hash_string;

    int f_count = 1; //从1开始，为了DSP优化
    int f_count_mod = f_count; //为了防止特征的数量超过MAX_FEATURE_NUMBER_IN_TILE
    int c_count = 0;
    int ff_count = 1;
    int ff_count_mod = ff_count;

    const enum RecognizeNoiseThreshold noise_threshold = context->noise_threshold;

#define LEFT_Y data[this_pos - 1]
#define TOP_Y data[top_line]
#define LEFT_TOP_Y data[top_line - 1]
#define TOP_RIGHT_Y data[top_line + 1]
#define THIS_Y data[this_pos]

    if (noise_threshold == kNoNoiseThreshold) {
        for (y = 1; y < 1079; y++) {
            for (x = 1; x < (frame->width - FEATURE_LENGTH); x++) {
                if (f_count >= MAX_FEATURE_NUMBER_IN_TILE)
                    break;
                this_pos = y * stride + x;
                top_line = this_pos - stride;
                //上及右上不同，则右边像素肯定也不是特征点
                if (TOP_Y ^ TOP_RIGHT_Y) {
                    x++;
                    continue;
                }
                if ((THIS_Y == LEFT_Y) ||
                    (LEFT_Y ^ LEFT_TOP_Y) ||
                    (LEFT_TOP_Y ^ TOP_Y))
                    continue;
                feature_item = feature_start + f_count_mod;
                feature_item->x = x;
                feature_item->y = y;
                hash_string = data + (feature_item->y) * stride + feature_item->x;
                //CRC CCITT-16
                //feature_item->hash_value = DoCrc(CRC_INIT, hash_string, FEATURE_LENGTH);
                feature_item->hash_value = calcHashKeyByStepXor13(hash_string, FEATURE_LENGTH);
                //此处使用的是临时缓冲区，内存不做拷贝，同步数据时再拷贝征点数据
                f_count++;
                f_count_mod = f_count % MAX_FEATURE_NUMBER_IN_TILE;
                x += FEATURE_LENGTH;
            }
            if (f_count >= MAX_FEATURE_NUMBER_IN_TILE)
                break;
        }
    }
    else{
        for (y = 1; y < 1079; y++) {
            for (x = 1; x < (frame->width - FEATURE_LENGTH); x++) {
                if (f_count >= MAX_FEATURE_NUMBER_IN_TILE){
                    break;
                }
                this_pos = y * stride + x;
                top_line = this_pos - stride;
                //利用条件提前退出做优化，不必做过多比较及运算

                //上及右上不同，则右边像素肯定也不是特征点
                if (!NOISE_SAME(TOP_Y ,TOP_RIGHT_Y,noise_threshold)) {
                    x++;
                    continue;
                }
                if ((NOISE_SAME(LEFT_Y, LEFT_TOP_Y, noise_threshold)) &&
                    (NOISE_SAME(LEFT_TOP_Y, TOP_Y, noise_threshold)) &&
                    (NOISE_SAME(LEFT_Y, TOP_Y, noise_threshold))&&
                    (NOISE_SAME(LEFT_Y, TOP_RIGHT_Y, noise_threshold))&&
                    (NOISE_SAME(LEFT_TOP_Y, TOP_RIGHT_Y, noise_threshold))&&
                    (! NOISE_SAME(LEFT_Y, THIS_Y, noise_threshold))){
                    feature_item = feature_start + f_count_mod;
                    feature_item->x = x;
                    feature_item->y = y;
                    hash_string = data + (feature_item->y) * stride + feature_item->x;
                    feature_item->hash_value = calcHashKeyByStepXor13(hash_string, FEATURE_LENGTH);
                    //此处使用的是临时缓冲区，内存不做拷贝，同步数据时再拷贝征点数据
                    f_count++;
                    f_count_mod = f_count % MAX_FEATURE_NUMBER_IN_TILE;
                    x += FEATURE_LENGTH;
                }
                if (f_count >= MAX_FEATURE_NUMBER_IN_TILE){
                    printf("break\n");
                    break;
                }

            }
        }
    }

    //防止f_count大于MAX_FEATURE_NUMBER_IN_TILE
    f_count = f_count > MAX_FEATURE_NUMBER_IN_TILE ? MAX_FEATURE_NUMBER_IN_TILE : f_count;
    FeatureIndexItem *index_item;
    FeatureItem * dst_item;
    for(int j=0;j<f_count;j++) {
        if (ff_count >= MAX_FEATURE_NUMBER_IN_TILE)
            break;
        //复制到目标表中
        feature_item = feature_start + j;
        hash_string = data + (feature_item->y) * stride + feature_item->x;
        dst_item = feature_table->feature_item + ff_count_mod;
        *dst_item = * feature_item;
        memcpy(dst_item->feature, hash_string, FEATURE_LENGTH);
        ff_count ++;
        ff_count_mod = ff_count % MAX_FEATURE_NUMBER_IN_TILE;

        index_item = index_table->feature_index_item + dst_item->hash_value;
        //记录特征的hash索引到index_table中
        if(0==index_item->feature_index){
            //没有冲突，记录索引
            index_item->feature_index = (unsigned short)(ff_count - 1);
        }else{
            c_count++;
            dst_item->collision_flag = 1; //标记该特征的HASH碰撞标记为1
            //(feature_table->feature_item + index_item->feature_index)->collision_flag = 1; //标记被碰撞的特征的HASH碰撞标记为1
        }
    }
    *feature_count = ff_count;
    *collision_count = c_count;
    return 0;
}

//发现当前帧和参考帧中匹配的所有特征
int findMatchedFeature(int feature_count_cur,FeatureTable *feature_table_cur, \
                       FeatureTable *feature_table_ref, \
                       FeatureIndexTable *index_table_ref, \
                       enum RecognizeNoiseThreshold noise_threshold, \
                       int *matched_count, MatchedFeatureTable *matched_table)
{
    int success = -1;
    int ii;
    int m_feature_count = 0;
    unsigned short hash_value;
    int crt_x;
    int crt_y;
    unsigned char *y_cur;
    int index;
    int ref_x;
    int ref_y;
    FeatureItem *feature_item_ref;
    FeatureIndexItem *index_item_ref;

    //因为matched_table相当与是这个函数使用的局部变量，这个函数要被调用多次，需要在这里清0，然后进行统计
    memset(matched_table, 0, sizeof(MatchedFeatureTable));
    //int feature_real = (feature_count_cur+1)<MAX_FEATURE_NUMBER_IN_TILE ? (feature_count_cur+1):MAX_FEATURE_NUMBER_IN_TILE;
    int feature_real = feature_count_cur;
    //对当前Tile中的每一个特征与参考Tile进行比对
    for (ii=1;ii<feature_real;ii++)
    {
        //从当前Tile的特征数组中得到一个特征
        FeatureItem *feature_item_cur = &(feature_table_cur->feature_item[ii]);
        if(1 == feature_item_cur->collision_flag){
            continue; //如果该特征与本帧中的其它特征发生了HASH碰撞，对该特征不做特征匹配计算
        }
        hash_value = feature_item_cur->hash_value;  //特征的Hash值
        crt_x = feature_item_cur->x; //特征的x坐标
        crt_y = feature_item_cur->y; //特征的y坐标
        //得到特征的亮度值数组，准备与参考Tile中的特征进行精确比对，避免Hash值冲突带来的影响
        y_cur = feature_item_cur->feature;

        //根据特征的Hash值可以迅速找到参考Tile中的匹配特征
        index_item_ref = &(index_table_ref->feature_index_item[hash_value]);
        //如果这个hash值对应的特征索引不为0，说明至少有一个参考Tile中存在的特征与这个特征的hash值一致，
        //有可能有一个是匹配特征，还需要进行精确验证
        //下面逐一进行精确比对
        if(0!=index_item_ref->feature_index){
            index = index_item_ref->feature_index;
            feature_item_ref = &(feature_table_ref->feature_item[index]);
            ref_x = feature_item_ref->x;
            ref_y = feature_item_ref->y;
            if ((ref_x == crt_x) && (ref_y == crt_y)) {
                //如果MV是(0,0)，直接跳过即可，因为我们有unchanged这种宏块的候选编码模式
                continue;
            } else{
                //如果MV不是(0,0)，需要进一步进行精确比对，排除hash冲突带来的影响
                unsigned char *y_ref = feature_item_ref->feature;

                // 噪音屏蔽范围为kNoNoiseThreshold时直接对比内存，否则逐像素对比
                if (noise_threshold == kNoNoiseThreshold) {
                    if(0 == memcmp(y_ref, y_cur, FEATURE_LENGTH)) {
                        //如果能够精确匹配，说明在当前Tile和参考Tile中存在一对匹配特征
                        //记录该特征的信息到列表中
                        MatchedFeatureItem *mf_item = &(matched_table->matched_feature_item[m_feature_count]);
                        mf_item->mv_x = ref_x - crt_x;
                        mf_item->mv_y = ref_y - crt_y;

#if RECORD_INFORMATION_FRAME_NUM
                        mf_item->local_x = crt_x;
                        mf_item->local_y = crt_y;
                        mf_item->hash_value = hash_value;
#endif

                        m_feature_count++;
                    }
                }
                else {
                    //精确对比特征的Y像素值
                    int is_same = 1;
                    for (int jj = 0; jj < FEATURE_LENGTH; jj++) {
                        if (!NOISE_SAME(*y_ref, *y_cur, noise_threshold)) {
                            is_same = 0;
                            break;
                        }
                        y_ref++;
                        y_cur++;
                    }

                    if (is_same == 1) {
                        //如果能够精确匹配，说明在当前Tile和参考Tile中存在一对匹配特征
                        //记录该特征的信息到列表中
                        MatchedFeatureItem *mf_item = &(matched_table->matched_feature_item[m_feature_count]);
                        mf_item->mv_x = ref_x - crt_x;
                        mf_item->mv_y = ref_y - crt_y;

#if RECORD_INFORMATION_FRAME_NUM
                        mf_item->local_x = crt_x;
                        mf_item->local_y = crt_y;
                        mf_item->hash_value = hash_value;
#endif

                        m_feature_count++;
                    }
                }
            }
        }
    }
    success = 0;
    *matched_count = m_feature_count;
    return success;
}


int findPossibleGlobalMv(int feature_count_cur, FeatureTable *feature_table_cur, \
                         FeatureTable *feature_table_ref, \
                         enum RecognizeNoiseThreshold noise_threshold, \
                         FeatureIndexTable *index_table_ref, \
                         int *matched_count, MatchedFeatureTable *matched_table, \
                         FoundMvTable *found_mv_table, struct GlobalMvMeta *global_mv)
{
    int success = -1;
    //利用参考Tile和当前Tile的特征信息找到匹配的特征集合
    int m_count = 0;
    FoundMvItem *found_mv_item;
    FoundMvItem *new_found_mv_item;
    MatchedFeatureItem *match_feature_item;
    int mv_x;
    int mv_y;

    success = findMatchedFeature(feature_count_cur, feature_table_cur, \
                                 feature_table_ref, index_table_ref, noise_threshold, \
                                 &m_count, matched_table);
    *matched_count = m_count;
    if(-1 == success)
    {
        printf("FindMatchedFeature has some problem.\n");
        return success;
    }
    //因为found_mv_table相当与是这个函数使用的局部变量，这个函数要被调用多次，需要在这里清0，然后进行统计
    memset(found_mv_table, 0, sizeof(FoundMvTable));

    int f_mv_count = 0;
    //逐个统计匹配特征对向量cv_matched_features_cur中各种mv出现的次数
    int ii,jj,kk;
    int found_same_mv_flag;
    for (jj=0;jj<m_count;jj++)
    {
        //统计各种MV出现的频次
        found_same_mv_flag = -1;
        match_feature_item = &(matched_table->matched_feature_item[jj]);
        mv_x = match_feature_item->mv_x;
        mv_y = match_feature_item->mv_y;
        //在found_mv_table中统计次数
        for (ii = 0; ii < f_mv_count; ii++) {
            found_mv_item = &(found_mv_table->found_mv_item[ii]);
            //统计相同的mv
            if ((mv_x == found_mv_item->mv_x) && (mv_y == found_mv_item->mv_y)) {
                found_mv_item->count ++;
                found_same_mv_flag = 1;
                break;
            }
        }
        //若果在mv统计向量cv_mv_meta中未发现，这是一个新出现的mv，插入到vector中
        if(-1 == found_same_mv_flag)
        {
            new_found_mv_item = &(found_mv_table->found_mv_item[f_mv_count]);
            new_found_mv_item->mv_x = mv_x;
            new_found_mv_item->mv_y = mv_y;
            new_found_mv_item->count = 1;
            f_mv_count++;
        }
    }

    //找到出现最频繁的mv，它有可能就是全局mv
    int count = 0;
    int index = -1;
    int valid_mv = 0;
    for(kk=0;kk<f_mv_count;kk++)
    {
        found_mv_item = &(found_mv_table->found_mv_item[kk]);
        if(count < found_mv_item->count)
        {
            count = found_mv_item->count;
            index = kk;
        }

        // 超过GLOBAL_MIX_MV_THRESHOLD_TO_MATCHED_FEATURES的MV，才做与最高频次MV的占比计算
        if (found_mv_item->count > GLOBAL_MIX_MV_THRESHOLD_TO_MATCHED_FEATURES) {
            valid_mv += found_mv_item->count;
        }

        //调试
        //printf("(%d, %d, %d)\n", found_mv_item->mv_x, found_mv_item->mv_y, found_mv_item->count);
    }
    //if((count > m_count*GLOBAL_MV_RATIO_THRESHOLD_TO_MATCHED_FEATURES)&&
    //   (count > feature_count_cur*GLOBAL_MV_RATIO_THRESHOLD_TO_ALL_FEATURES))
    if((count > valid_mv*GLOBAL_MV_RATIO_THRESHOLD_TO_MATCHED_FEATURES)&&(count > GLOBAL_MV_THRESHOLD_TO_MATCHED_FEATURES_NUM))
    {
        //如果出现频次最高的mv占比超过了门限值，认为它是全局mv
        found_mv_item = &(found_mv_table->found_mv_item[index]);
        global_mv->valid_flag = 1;
        global_mv->x = found_mv_item->mv_x;
        global_mv->y = found_mv_item->mv_y;
        //printf("total feature count:%d, matched feature count:%d, max matched feature count:%d\n", feature_count_cur, m_count, count);
    }
    else
    {
        //否则，认为是噪声mv，该Tile与参考Tile不存在全局mv
        global_mv->valid_flag = -1;
        global_mv->x = 0;
        global_mv->y = 0;
        //VLOGD("No possible global MV..\n");
    }

    success = 0;
    return success;
}


int findFeatureInfoAndGlobalMv(struct RecognizeContext* const context){
    int success = -1;
    int find_feature_index_table_flag = -1;
    int *ref_frame_id = &(context->frame_id_of_one_frame); //先标记参考帧用frame_id_of_one_frame
    FeatureTable *feature_table_ref = &(context->feature_meta_in_one_frame->feature_table);
    FeatureIndexTable *index_table_ref = &(context->feature_meta_in_one_frame->index_table);
    int *cur_frame_id = &(context->frame_id_of_another_frame); //先标记当前帧用frame_id_of_another_frame
    FeatureTable *feature_table_cur = &(context->feature_meta_in_another_frame->feature_table);
    FeatureIndexTable *index_table_cur = &(context->feature_meta_in_another_frame->index_table);
    if(context->frame_id_of_one_frame == context->refer_frame->frame_id){
        //这种情况下，符合缺省的情况
        //在another中存在的特征索引表不需要了，当前帧的特征索引表数据可以存储在其中
        //找到了参考帧的特征索引表，说明参考帧已经做过了特征搜索了
        find_feature_index_table_flag = 1;
    }else if(context->frame_id_of_another_frame == context->refer_frame->frame_id){
        //找到了参考帧的特征索引表，说明参考帧已经做过了特征搜索了
        find_feature_index_table_flag = 1;
        ref_frame_id = &(context->frame_id_of_another_frame);
        feature_table_ref = &(context->feature_meta_in_another_frame->feature_table);
        index_table_ref = &(context->feature_meta_in_another_frame->index_table);
        //在one中存在的特征索引表不需要了，当前帧的特征索引表数据可以存储在其中
        cur_frame_id = &(context->frame_id_of_one_frame);
        feature_table_cur = &(context->feature_meta_in_one_frame->feature_table);
        index_table_cur = &(context->feature_meta_in_one_frame->index_table);
    }
    if(-1 == find_feature_index_table_flag){
        //参考帧没有进行过特征搜索，对参考帧进行特征搜索
        int feature_count_ref = -1;
        int c_count_ref = 0;
        //参考帧特征索引表需要清零
        memset(index_table_ref, 0, sizeof(FeatureIndexTable));
        memset(feature_table_ref, 0, sizeof(FeatureTable));
        success = findAllFeaturesInFrame(context,context->refer_frame, &feature_count_ref, &c_count_ref, feature_table_ref,
                                         index_table_ref);
        *ref_frame_id = context->refer_frame->frame_id;
    }
    //对当前帧进行特征搜索
    int feature_count_cur = -1;
    int c_count_cur = 0;
    //当前帧特征索引表需要清零
    memset(index_table_cur, 0, sizeof(FeatureIndexTable));
    memset(feature_table_cur, 0, sizeof(FeatureTable));
    success = findAllFeaturesInFrame(context,context->current_frame, &feature_count_cur, &c_count_cur, feature_table_cur,
                                     index_table_cur);
    *cur_frame_id = context->current_frame->frame_id;
    //3.2.4 进行特征匹配，查找MV
    //分配临时变量，准备做特征匹配
    int matched_count = -1;
    MatchedFeatureTable matched_table = {0};
    FoundMvTable found_mv_table = {0};
    //继续做ME等识别处理
    //找到匹配的特征对
    //估计当前Tile相对于参考Tile的全局mv
    struct GlobalMvMeta global_mv;
    success = findPossibleGlobalMv(feature_count_cur, feature_table_cur, \
                                   feature_table_ref, context->noise_threshold, \
                                   index_table_ref, \
                                   &matched_count, &matched_table, &found_mv_table, &global_mv);
    if (1 == global_mv.valid_flag) {
        printf("The global mv is (%d, %d)\n", global_mv.x, global_mv.y);
        context->global_mv.x = global_mv.x;
        context->global_mv.y = global_mv.y;
        context->global_mv.valid_flag = 1; //设置全局MV有效
    }
    else {
        //printf("The global mv not exsit.\n");
        context->global_mv.x = 0;
        context->global_mv.y = 0;
        context->global_mv.valid_flag = -1; //设置全局MV有效
    }

#if RECORD_INFORMATION_FRAME_NUM
    if (context->current_frame->frame_id == RECORD_INFORMATION_FRAME_NUM && matched_count > 0) {
        FILE* allFeatureRecoder = fopen("./AllFeatures.txt", "wb");
        if (allFeatureRecoder) {
            for (int i = 0; i < MAX_FEATURE_NUMBER_IN_TILE; i++) {
                if (feature_table_cur->feature_item[i].collision_flag == 0) {
                    fprintf(allFeatureRecoder, "%d, %d, %d, %02x\n", i,
                            feature_table_cur->feature_item[i].x,
                            feature_table_cur->feature_item[i].y,
                            feature_table_cur->feature_item[i].hash_value);
                }
            }
            fclose(allFeatureRecoder);
        }

        FILE* matchedFeatureRecoder = fopen("./MatchedFeatures.txt", "wb");
        if (matchedFeatureRecoder) {
            for (int i = 0; i < matched_count; i++) {
                fprintf(matchedFeatureRecoder, "%d, %d, %d, %d, %d, %02x\n", i,
                        matched_table.matched_feature_item[i].local_x, \
						matched_table.matched_feature_item[i].local_y, \
						matched_table.matched_feature_item[i].mv_x, \
						matched_table.matched_feature_item[i].mv_y, \
						matched_table.matched_feature_item[i].hash_value);
            }
            fclose(matchedFeatureRecoder);
        }
    }
#endif

    //MV匹配过程完成，
    return success;
}

//比较当前帧和参考帧运动补偿后（位移全局MV后）的宏块是否相同
int fastCompareMBWithGlobalMVInReference(int mb_x, int mb_y, int mb_width, int mb_height, \
                                     int mv_x, int mv_y, \
                                     struct Frame *cur, struct Frame *ref, \
                                     enum RecognizeNoiseThreshold noise_threshold) {
    DEF_COLOR_FORMAT_LOCAL_VAR(cur->stride);
    //在合理范围之内，不必判断如果使用全局MV做补偿是否在参考Tile的范围之外
    //TODO(ssurui)可优化，减少计算次数
    unsigned char *p_data1_y = MAP_MB_DATA(1,cur->data_component1,mb_y,mb_x);
    unsigned char *p_data1_u = MAP_MB_DATA(2,cur->data_component2,mb_y,mb_x);
    unsigned char *p_data1_v = MAP_MB_DATA(3,cur->data_component3,mb_y,mb_x);
    unsigned char *p_data2_y = MAP_MB_DATA(1,ref->data_component1,mb_y,mb_x) + mv_y  * local_data1_stride + mv_x;
    unsigned char *p_data2_u = MAP_MB_DATA(2,ref->data_component2,mb_y,mb_x);
    unsigned char *p_data2_v = MAP_MB_DATA(3,ref->data_component3,mb_y,mb_x);
    if (IS_YUV420P){
        p_data2_u += mv_y / 2 * local_data2_stride + mv_x / 2;
        p_data2_v += mv_y / 2 * local_data2_stride + mv_x / 2;
    }else{
        p_data2_u += mv_y * local_data2_stride + mv_x;
        p_data2_v += mv_y * local_data2_stride + mv_x;
    }

    if (kNoNoiseThreshold == noise_threshold) {
        int y;
        for (y = 0; y < MB_HEIGHT; y ++){
            if (0 != memcmp(p_data1_y, p_data2_y, local_data1_MB_WIDTH))
                return -1;
            p_data1_y += local_data1_stride;
            p_data2_y += local_data1_stride;
            if (IS_YUV420P && (y % 2 != 0))
                continue;
            if ((0 != memcmp(p_data1_u, p_data2_u, local_data2_MB_WIDTH)) ||
                (0 != memcmp(p_data1_v, p_data2_v, local_data3_MB_WIDTH))) {
                return -1;
            }
            p_data1_u += local_data2_stride;
            p_data2_u += local_data2_stride;
            p_data1_v += local_data3_stride;
            p_data2_v += local_data3_stride;
        }
    }
    else {
        unsigned char *p_mb_y_1 = p_data1_y;
        unsigned char *p_mb_u_1 = p_data1_u;
        unsigned char *p_mb_v_1 = p_data1_v;
        unsigned char *p_mb_y_2 = p_data2_y;
        unsigned char *p_mb_u_2 = p_data2_u;
        unsigned char *p_mb_v_2 = p_data2_v;
        int local_data1_jump = local_data1_stride - local_data1_MB_WIDTH;
        int local_data2_jump = local_data2_stride - local_data2_MB_WIDTH;
        int local_data3_jump = local_data3_stride - local_data3_MB_WIDTH;

        for (int j = 0; j < MB_HEIGHT; j++) {
            for (int i = 0; i < MB_WIDTH; i++) {
                if (!NOISE_SAME(*p_mb_y_1, *p_mb_y_2, noise_threshold))
                    return -1;
                p_mb_y_1++;
                p_mb_y_2++;
                if (IS_YUV420P && ((j %2 != 0)||(i %2 != 0)))
                    continue;
                if (!NOISE_SAME(*p_mb_u_1, *p_mb_u_2, noise_threshold) ||
                    !NOISE_SAME(*p_mb_v_1, *p_mb_v_2, noise_threshold)) {
                    return -1;
                }

                p_mb_u_1++;
                p_mb_v_1++;

                p_mb_u_2++;
                p_mb_v_2++;
            }
            p_mb_y_1 += local_data1_jump;
            p_mb_y_2 += local_data1_jump;

            if (!IS_YUV420P || (j % 2 == 0)){
                p_mb_u_1 += local_data2_jump;
                p_mb_u_2 += local_data2_jump;

                p_mb_v_1 += local_data3_jump;
                p_mb_v_2 += local_data3_jump;
            }
        }
    }
    return 1;
}


int markMvMatchedMb(struct RecognizeContext *const context) {
    if(1 != context->global_mv.valid_flag){
        //运动向量不存在，直接退出
        return 0;
    }
    //printf("frame_id: %d has a global mv: (%d, %d)\n", context->current_frame->frame_id, context->global_mv.x, context->global_mv.y);
//    int mb_in_column = context->mb_in_column;
//    int mb_in_row = context->mb_in_row;
    int mb_stride = context->mb_stride;
    int mv_x = context->global_mv.x;
    int mv_y = context->global_mv.y;

    //由于全局运动的关系，有可能的宏块位于当前帧和参考帧的重叠范围
    //先求出这个重叠范围
    int x_start = (0 > -context->global_mv.x) ? 0 : -context->global_mv.x;
    int x_end = (context->frame_width < context->frame_width - context->global_mv.x) ? context->frame_width : context->frame_width - context->global_mv.x;
    int y_start = (0 > -context->global_mv.y) ? 0 : -context->global_mv.y;
    int y_end = (context->frame_height < context->frame_height - context->global_mv.y) ? context->frame_height : context->frame_height - context->global_mv.y;
    //将范围转化为宏块范围
    int mb_x_start = (x_start == 0) ? 0 : (x_start-1)/MB_WIDTH + 1;
    int mb_x_end = (x_end-1)/MB_WIDTH ;
    int mb_y_start = (y_start == 0) ? 0 : (y_start-1)/MB_HEIGHT + 1;
    int mb_y_end = (y_end-1)/MB_HEIGHT;

    //对属于可能的运动补偿宏块进行重新标记
    int same_flag;
    enum MBType *mb_type = NULL;
    for(int jj = mb_y_start; jj < mb_y_end;jj++) {
        for (int ii = mb_x_start; ii < mb_x_end; ii++) {
            //逐个宏块检查
            mb_type = context->mbtype + jj*mb_stride + ii;
            if((DEFAULT == *mb_type)&&(TEXT != *mb_type)){
                //只有这两部分宏块才有可能是运动补偿宏块
                same_flag = fastCompareMBWithGlobalMVInReference(ii, jj, MB_WIDTH, MB_HEIGHT, mv_x, mv_y, context->current_frame, context->refer_frame, context->noise_threshold);
                if(1 == same_flag) {
                    *mb_type = MVMATCHED;
                    context->mvmatched_count ++;
                }
            }
        }
    }
    return 0;
}

///检测当前帧相对于参考帧的运动信息
int generateMotionDetect(struct Frame * current_frame, struct Frame *refer_frame, struct RecognizeContext *const context){
    int success = -1;
    //第一步，简单赋值
    context->current_frame = (struct Frame *)current_frame;
    context->refer_frame = (struct Frame *)refer_frame;  //赋值参考帧指针
    context->global_mv.valid_flag = -1; //初始化全局运动向量为无效
    context->global_mv.x = 0;
    context->global_mv.y = 0;
    context->mvmatched_count = 0;
    if(NULL != context->refer_frame){
        success = findFeatureInfoAndGlobalMv(context);
        success = markMvMatchedMb(context);
    }

    return success;
}


/*
int nbcGenerateSmart_Ex(unsigned char *p_data_y, int stride, int mb_height, int mb_width,int max) {
    int nbc_block = 0;  //主颜色个数
    //int value_bak[3] = {0,0,0}; //查找备份
    int color;  //临时变量
    const int nbc_diff = NBC_DIFF;  //如果一种颜色与主颜色值的差小于该值，认为这种颜色也是主颜色
    int m, n, i, j;  //临时循环变量
    const int SIZE = 256;
    int same_value_count[257] = {0};
    unsigned char *p_data_tmp = p_data_y;
    int stride_diff = stride - mb_width;

    //统计Y分量的颜色直方图，结果记录在value数组中
    for (m = 0; m < mb_height; m++) {
        for (n = 0; n < mb_width; n++) {
            color = *p_data_tmp;
            same_value_count[color]++;
            ++p_data_tmp;
        }
        p_data_tmp += stride_diff;
    }

    //填充颜色范围是0-255，下标256初始化成0
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
            //从前到后遍历找插入位置
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
    //返回主颜色所占的像素个数
    return nbc_block;
}
*/

/*
int nhgGenerateSmart_Ex(unsigned char *p_data_y, int stride, int mb_height, int mb_width, int max_value) {
    int i, j;
    int nhg_block = 0; //高梯度值像素的个数
    int grads = 0; //梯度，临时变量
    unsigned char *p_data_tmp = p_data_y;
    int stride_diff = stride - mb_width;
    int mb_width_last_p = mb_width - 1;

    //梯度值计算，只要该像素与周围8个像素的Y分量的差值有一个大于nhg_diff，该像素就属于高梯度值像素
    for (i = 0; i < mb_height; ++i) {
        for (j = 0; j < mb_width; ++j,++p_data_tmp) {
            grads = abs(*p_data_tmp - *(p_data_tmp +1));
            if (grads > NGH_DIFF) {
                if(j != mb_width_last_p){
                    //非最后一列，下一个点可以不用检查直接满足条件
                    ++nhg_block;
                    ++j;
                    ++p_data_tmp;
                }
                //if (++nhg_block　ｖ
                //
                // > max_value).

                //    return nhg_block;
                ++nhg_block;
                continue;
            }
            grads = abs(*p_data_tmp - *(p_data_tmp + stride));
            if (grads > NGH_DIFF) {
                // if (++nhg_block > max_value)
                //  return nhg_block;
                ++nhg_block;
                continue;
            }
            grads = abs(*p_data_tmp - *(p_data_tmp - stride));
            if (grads > NGH_DIFF) {
                //if (++nhg_block > max_value)
                //return nhg_block;
                ++nhg_block;
                continue;
            }
            grads = abs(*p_data_tmp - *(p_data_tmp - 1));
            if (grads > NGH_DIFF) {
                //if (++nhg_block > max_value)
                //return nhg_block;
                ++nhg_block;
                continue;
            }
        }
        p_data_tmp += stride_diff;
    }
    return nhg_block;
}

//检测文字块
bool decideTextBlock(unsigned char *p_data_y, int local_data1_stride){
    int mb_nhg_count = 0;
    int mb_nbc_count = 0;
    bool ret_val = false;
    mb_nbc_count = nbcGenerateSmart_Ex(p_data_y, local_data1_stride, MB_HEIGHT, MB_WIDTH, NBC_TH0);
    mb_nhg_count = nhgGenerateSmart_Ex(p_data_y, local_data1_stride, MB_HEIGHT, MB_WIDTH, NGH_TH);
    if (mb_nhg_count < NGH_TH){
        ret_val = false;
    } else{
        if (mb_nbc_count > 160) {
            ret_val = true;  //kText
            if (mb_nhg_count < 20) {
                ret_val = false;
            }
        }else {
            ret_val = true;  //kText
            if (mb_nhg_count < 17) {
                ret_val = false;
            }
        }
    }
    return ret_val;
}
*/
/*
void detectTextMB(struct Frame *cur_frm, struct RecognizeContext *const context){
    DEF_COLOR_FORMAT_LOCAL_VAR(cur_frm->stride);
    context->text_count = 0;
    int mb_nhg_count = 0;
    int mb_nbc_count = 0;
    enum MBType *mb_type = NULL;
    for (int i = 0; i < context->mb_in_column; ++i) {
        for (int j = 0; j < context->mb_in_row; ++j) {
            mb_type = context->mbtype + i*context->mb_stride + j;
            if (DEFAULT == *mb_type){
                unsigned char *p_data_y = MAP_MB_DATA(1,cur_frm->data_component1,i,j);//宏块左上角像素的Y分量的数据指针
                mb_nbc_count = nbcGenerateSmart_Ex(p_data_y, local_data1_stride, local_data1_MB_HEIGHT, local_data1_MB_WIDTH, NBC_TH0);
                mb_nhg_count = nhgGenerateSmart_Ex(p_data_y, local_data1_stride, local_data1_MB_HEIGHT, local_data1_MB_WIDTH, NGH_TH);

                enum MBType block_mode_left = *(mb_type - 1);
                enum MBType block_mode_up = *(mb_type - context->mb_stride);
                enum MBType block_mode_upleft = *(mb_type - context->mb_stride - 1);

                if (mb_nhg_count <= NGH_TH) {
                    continue;
                }
                else{
                    if((block_mode_left == TEXT) && (block_mode_up == TEXT) && (block_mode_upleft == TEXT)){
                        if (mb_nbc_count > NBC_TH1){
                            *mb_type = TEXT;
                            context->text_count ++;
                        }
                    }else{
                        if (mb_nbc_count > NBC_TH0){
                            *mb_type = TEXT;
                            context->text_count ++;
                        }
                    }
                }
            }
        }
    }
}
*/


bool getSamePixelRatio(unsigned char *p_data1_y, unsigned char *p_data1_u, unsigned char *p_data1_v, \
                        unsigned char *p_data2_y, unsigned char *p_data2_u, unsigned char *p_data2_v, \
                        struct Frame *frm){
    DEF_COLOR_FORMAT_LOCAL_VAR(frm->stride);
    int y_same_num = 0;
    int u_same_num = 0;
    int v_same_num = 0;
    int total_pixels = 0;
    int stride1 = local_data1_stride;
    int stride2 = local_data2_stride;
    int stride3 = local_data3_stride;
    //比较分量
    if (IS_YUV420P){
        total_pixels = local_data1_MB_PIXEL_COUNT + local_data2_MB_PIXEL_COUNT + local_data3_MB_PIXEL_COUNT;
        for (int y = 0; y < MB_HEIGHT; y++) {
            for (int x = 0; x < MB_WIDTH; x++) {
                if (abs(*(p_data1_y + x) - *(p_data2_y + x)) <= 1){
                    y_same_num ++;
                }
                if((y % 2 != 0) || (x % 2 != 0)){
                    continue;
                }
                if (abs(*(p_data1_u + x/2) - *(p_data2_u + x/2)) <= 1){
                    u_same_num ++;
                }
                if (abs(*(p_data1_v + x/2) - *(p_data2_v + x/2)) <= 1){
                    v_same_num ++;
                }
            }
            p_data1_y += stride1;
            p_data2_y += stride1;
            if (y % 2 == 0){
                p_data1_u += stride2;
                p_data2_u += stride2;
                p_data1_v += stride3;
                p_data2_v += stride3;
            }
        }
    } else{
        total_pixels = MB_HEIGHT*MB_WIDTH*3;
        for (int y = 0; y < MB_HEIGHT; y++) {
            for (int x = 0; x < MB_WIDTH; x++) {
                if (abs(*(p_data1_y + x) - *(p_data2_y + x)) <= 1){
                    y_same_num ++;
                }
                if (abs(*(p_data1_u + x) - *(p_data2_u + x)) <= 1){
                    u_same_num ++;
                }
                if (abs(*(p_data1_v + x) - *(p_data2_v + x)) <= 1){
                    v_same_num ++;
                }
            }
            p_data1_y += stride1;
            p_data2_y += stride1;
            p_data1_u += stride2;
            p_data2_u += stride2;
            p_data1_v += stride3;
            p_data2_v += stride3;
        }
    }
    return (y_same_num + u_same_num + v_same_num == total_pixels);
}

//检测IBC
void detectIntraCopyBlock(struct Frame *cur_frm, struct RecognizeContext *const context){
    DEF_COLOR_FORMAT_LOCAL_VAR(cur_frm->stride);

    context->lbc_count = 0;
    context->ubc_count = 0;
    int ref_mb_x = 0;
    int ref_mb_y = 0;
    bool up_same_flag = 0;
    bool left_same_flag = 0;
    enum MBType *mb_type = NULL;
    for (int i = 0; i < context->mb_in_column; ++i) {
        for (int j = 0; j < context->mb_in_row; ++j) {
			if (i == 1 && j == 0)
				ref_mb_x = 0;

            ref_mb_x = j - 1;
            ref_mb_y = i - 1;
            mb_type = context->mbtype + i * context->mb_stride + j;
            unsigned char *p_data_y = MAP_MB_DATA(1,cur_frm->data_component1,i,j);//宏块左上角像素的Y分量的数据指针
            unsigned char *p_data_u = MAP_MB_DATA(2,cur_frm->data_component2,i,j); //宏块左上角像素的U分量的数据指针
            unsigned char *p_data_v = MAP_MB_DATA(3,cur_frm->data_component3,i,j); //宏块左上角像素的V分量的数据指针
            unsigned char *p_up_data_y = MAP_MB_DATA(1,cur_frm->data_component1,ref_mb_y,j);//宏块左上角像素的Y分量的数据指针
            unsigned char *p_up_data_u = MAP_MB_DATA(2,cur_frm->data_component2,ref_mb_y,j); //宏块左上角像素的U分量的数据指针
            unsigned char *p_up_data_v = MAP_MB_DATA(3,cur_frm->data_component3,ref_mb_y,j); //宏块左上角像素的V分量的数据指针
            unsigned char *p_left_data_y = MAP_MB_DATA(1,cur_frm->data_component1,i,ref_mb_x);//宏块左上角像素的Y分量的数据指针
            unsigned char *p_left_data_u = MAP_MB_DATA(2,cur_frm->data_component2,i,ref_mb_x); //宏块左上角像素的U分量的数据指针
            unsigned char *p_left_data_v = MAP_MB_DATA(3,cur_frm->data_component3,i,ref_mb_x); //宏块左上角像素的V分量的数据指针

            up_same_flag = getSamePixelRatio(p_data_y,p_data_u,p_data_v,p_up_data_y,p_up_data_u,p_up_data_v, cur_frm);
            left_same_flag = getSamePixelRatio(p_data_y,p_data_u,p_data_v,p_left_data_y,p_left_data_u,p_left_data_v, cur_frm);

            if (left_same_flag){
                *mb_type = LBC;
                context->lbc_count ++;
            } else if (up_same_flag){
                *mb_type = UBC;
                context->ubc_count ++;
            }
        }
    }
}


