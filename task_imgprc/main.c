#define _POSIX_C_SOURCE 200809L

#include <opencv2/core/core_c.h>
#include <opencv2/highgui/highgui_c.h>
#include <opencv2/imgproc/imgproc_c.h>

#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#define MAX_REGIONS             128
#define MAX_TRACKS              64

#define CV_RGB(r, g, b)         cvScalar((b), (g), (r), 0)

#define DWT_LEVELS              1
#define INIT_DWT_THRESH_INT     50
#define INIT_MIN_CONF_INT       40

#define WIN_W                   800
#define WIN_H                   600
#define POSITION_TICKS          1000

#define SEEK_STEP_US            (5LL * 1000000)
#define SEEK_BIG_FACTOR         6
#define FRAME_SKIP_MOD          2

#define REGION_MIN_WIDTH_PCT    0.05f
#define REGION_MIN_WIDTH_FLOOR  80
#define REGION_MAX_WIDTH_PCT    0.98f
#define REGION_MIN_HEIGHT_PCT   0.02f
#define REGION_MIN_HEIGHT_FLOOR 14
#define REGION_MAX_HEIGHT_PCT   0.12f
#define REGION_MAX_HEIGHT_FLOOR 40
#define REGION_MIN_ASPECT       3.0f
#define REGION_MAX_ASPECT       80.0f
#define REGION_TOP_LIMIT_PCT    0.65f
#define REGION_BOTTOM_LIMIT_PCT 0.995f
#define REGION_BOTTOM_BAND_PCT  0.78f
#define REGION_CENTER_TOL_PCT   0.5f
#define REGION_MIN_FULLNESS     0.25f
#define REGION_MIN_CONTOUR_LEN  8

#define BRIGHT_THRESHOLD        180
#define BRIGHT_RATIO_MIN        0.02f
#define BRIGHT_RATIO_MAX        0.98f
#define CONTRAST_MIN            0.35f

#define NMS_IOU_THRESH          0.4f

#define TRACK_IOU_MIN           0.3f
#define TRACK_HIT_EMIT          2
#define TRACK_HIT_STABLE        4
#define TRACK_MAX_MISSES        8
#define TRACK_EMA_KEEP          0.6f
#define TRACK_CONF_EMA_KEEP     0.7f

#define THRESH_BLOCK_SIZE       25
#define THRESH_C                8.0f
#define THRESH_HH_FACTOR        0.85f

#define MORPH_LH_TEXT_W         5
#define MORPH_LH_TEXT_H         3
#define MORPH_LH_CLEAN_W        3
#define MORPH_LH_CLEAN_H        3
#define MORPH_HL_TEXT_W         3
#define MORPH_HL_TEXT_H         7
#define MORPH_HL_CLEAN_W        2
#define MORPH_HL_CLEAN_H        2
#define MORPH_HH_TEXT_W         3
#define MORPH_HH_TEXT_H         3
#define MORPH_HH_CLEAN_W        2
#define MORPH_HH_CLEAN_H        2

#define LINE_KERNEL_W           51
#define LINE_KERNEL_H           3
#define LINE_DILATE_ITERS       2
#define CLEAN_KERNEL_W          2
#define CLEAN_KERNEL_H          2

struct TextRegion {
	CvRect rect;
	float  confidence;
	int    is_stable;
};

struct TextRegionList {
	struct TextRegion regions[MAX_REGIONS];
	int               count;
};

struct Track {
	CvRect rect;
	float  confidence;
	int    hits;
	int    misses;
	int    active;
	int    captured;
};

struct Tracker {
	struct Track tracks[MAX_TRACKS];
	int          count;
};

struct Video {
	AVFormatContext   *fmt;
	AVCodecContext    *dec;
	struct SwsContext *sws;
	AVFrame           *frame;
	AVFrame           *bgr;
	AVPacket          *pkt;
	uint8_t           *bgr_buf;
	int                video_stream;
	int                width;
	int                height;
	double             fps;
	AVRational         stream_tb;
	int64_t            duration_us;
	int64_t            cur_pts_us;
};

struct DwtLevel {
	CvMat *low;
	CvMat *high;
	CvMat *LL;
	CvMat *LH;
	CvMat *HL;
	CvMat *HH;
};

struct ThreshScratch {
	CvMat    *abs_detail;
	IplImage *normalized;
	IplImage *global_binary;
};

struct MorphSet {
	IplImage      *temp;
	IplConvKernel *k_text;
	IplConvKernel *k_clean;
};

struct ProcCtx {
	int                  W, H;
	int                  dW, dH;
	int                  levels;

	IplImage            *gray;
	CvMat               *gray32;

	struct DwtLevel     *dwt;

	struct ThreshScratch t_scratch;
	IplImage            *binLH, *binHL, *binHH;

	struct MorphSet      morphLH_set, morphHL_set, morphHH_set;
	IplImage            *morphLH, *morphHL, *morphHH;

	IplImage            *edge_map;
	IplImage            *edge_copy;
	IplConvKernel       *line_kernel;
	IplConvKernel       *clean_kernel;
	CvMemStorage        *storage;

	IplImage            *display;
	IplImage            *enh_gray;
	IplImage            *enh_binary;

	IplImage            *vchan;
	IplImage            *salience;

	unsigned char       *ppm_buf;
	size_t               ppm_buf_sz;
};

static double
now_seconds(void)
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec / 1e9;
}

static void
ensure_dir(const char *path)
{
	struct stat st;

	if (stat(path, &st) != 0)
		mkdir(path, 0755);
}

static void
safe_release_mat(CvMat **mat)
{
	if (mat && *mat) {
		cvReleaseMat(mat);
		*mat = NULL;
	}
}

static void
safe_release_image(IplImage **img)
{
	if (img && *img) {
		cvReleaseImage(img);
		*img = NULL;
	}
}

static void
safe_release_kernel(IplConvKernel **k)
{
	if (k && *k) {
		cvReleaseStructuringElement(k);
		*k = NULL;
	}
}

static float
rect_iou(CvRect a, CvRect b)
{
	int   ix1, iy1, ix2, iy2;
	int   ax2, ay2, bx2, by2;
	int   iw, ih;
	float inter, uni;

	ix1 = a.x > b.x ? a.x : b.x;
	iy1 = a.y > b.y ? a.y : b.y;
	ax2 = a.x + a.width;
	ay2 = a.y + a.height;
	bx2 = b.x + b.width;
	by2 = b.y + b.height;
	ix2 = ax2 < bx2 ? ax2 : bx2;
	iy2 = ay2 < by2 ? ay2 : by2;
	iw  = ix2 - ix1;
	ih  = iy2 - iy1;

	if (iw <= 0 || ih <= 0)
		return 0.0f;

	inter = (float)(iw * ih);
	uni   = (float)(a.width * a.height + b.width * b.height) - inter;

	return uni > 0.0f ? inter / uni : 0.0f;
}

static void
emit_frame(struct Video *v, IplImage *out)
{
	int y;

	sws_scale(v->sws, (const uint8_t *const *)v->frame->data, v->frame->linesize,
	          0, v->height, v->bgr->data, v->bgr->linesize);

	for (y = 0; y < v->height; ++y) {
		memcpy(out->imageData + y * out->widthStep,
		       v->bgr->data[0] + y * v->bgr->linesize[0],
		       v->width * 3);
	}

	if (v->frame->best_effort_timestamp != AV_NOPTS_VALUE) {
		v->cur_pts_us = av_rescale_q(v->frame->best_effort_timestamp,
		                             v->stream_tb, AV_TIME_BASE_Q);
	}
}

static int
video_open(struct Video *v, const char *path)
{
	unsigned       i;
	AVStream      *st;
	const AVCodec *codec;
	int            bytes;

	memset(v, 0, sizeof(*v));
	v->video_stream = -1;

	if (avformat_open_input(&v->fmt, path, NULL, NULL) < 0)
		return 0;
	if (avformat_find_stream_info(v->fmt, NULL) < 0)
		return 0;

	for (i = 0; i < v->fmt->nb_streams; ++i) {
		if (v->fmt->streams[i]->codecpar->codec_type == AVMEDIA_TYPE_VIDEO) {
			v->video_stream = (int)i;
			break;
		}
	}
	if (v->video_stream < 0)
		return 0;

	st    = v->fmt->streams[v->video_stream];
	codec = avcodec_find_decoder(st->codecpar->codec_id);
	if (!codec)
		return 0;

	v->dec = avcodec_alloc_context3(codec);
	if (!v->dec)
		return 0;
	if (avcodec_parameters_to_context(v->dec, st->codecpar) < 0)
		return 0;
	if (avcodec_open2(v->dec, codec, NULL) < 0)
		return 0;

	v->width  = v->dec->width;
	v->height = v->dec->height;
	v->fps    = av_q2d(st->avg_frame_rate);
	if (v->fps <= 0.0)
		v->fps = av_q2d(st->r_frame_rate);
	v->stream_tb   = st->time_base;
	v->duration_us = v->fmt->duration > 0 ? v->fmt->duration : 0;
	v->cur_pts_us  = 0;

	v->frame       = av_frame_alloc();
	v->bgr         = av_frame_alloc();
	v->pkt         = av_packet_alloc();
	if (!v->frame || !v->bgr || !v->pkt)
		return 0;

	bytes      = av_image_get_buffer_size(AV_PIX_FMT_BGR24, v->width, v->height, 1);
	v->bgr_buf = (uint8_t *)av_malloc(bytes);
	av_image_fill_arrays(v->bgr->data, v->bgr->linesize, v->bgr_buf,
	                     AV_PIX_FMT_BGR24, v->width, v->height, 1);

	v->sws = sws_getContext(v->width, v->height, v->dec->pix_fmt,
	                        v->width, v->height, AV_PIX_FMT_BGR24,
	                        SWS_BILINEAR, NULL, NULL, NULL);
	if (!v->sws)
		return 0;

	return 1;
}

static int
video_read(struct Video *v, IplImage *out)
{
	int ret;

	while (1) {
		ret = avcodec_receive_frame(v->dec, v->frame);
		if (ret == 0) {
			emit_frame(v, out);
			return 1;
		}
		if (ret != AVERROR(EAGAIN))
			return 0;

		ret = av_read_frame(v->fmt, v->pkt);
		if (ret < 0) {
			avcodec_send_packet(v->dec, NULL);
			ret = avcodec_receive_frame(v->dec, v->frame);
			if (ret < 0)
				return 0;
			emit_frame(v, out);
			return 1;
		}
		if (v->pkt->stream_index == v->video_stream)
			avcodec_send_packet(v->dec, v->pkt);
		av_packet_unref(v->pkt);
	}
}

static int
video_seek(struct Video *v, int64_t target_us)
{
	int64_t target_ts;

	if (target_us < 0)
		target_us = 0;
	if (v->duration_us > 0 && target_us > v->duration_us)
		target_us = v->duration_us;

	target_ts = av_rescale_q(target_us, AV_TIME_BASE_Q, v->stream_tb);
	if (av_seek_frame(v->fmt, v->video_stream, target_ts, AVSEEK_FLAG_BACKWARD) < 0)
		return 0;

	avcodec_flush_buffers(v->dec);
	v->cur_pts_us = target_us;
	return 1;
}

static void
seek_reset(struct Video *v, struct Tracker *tk, int64_t target_us)
{
	video_seek(v, target_us);
	memset(tk, 0, sizeof(*tk));
}

static void
video_close(struct Video *v)
{
	if (v->sws)
		sws_freeContext(v->sws);
	if (v->bgr_buf)
		av_free(v->bgr_buf);
	if (v->frame)
		av_frame_free(&v->frame);
	if (v->bgr)
		av_frame_free(&v->bgr);
	if (v->pkt)
		av_packet_free(&v->pkt);
	if (v->dec)
		avcodec_free_context(&v->dec);
	if (v->fmt)
		avformat_close_input(&v->fmt);
}

static void
dwt_row_pass(const CvMat *in, CvMat *low, CvMat *high)
{
	const float inv_sqrt2 = 1.0f / sqrtf(2.0f);
	int         r         = in->rows;
	int         c         = in->cols;
	int         i, j;
	float       a, b;
	float      *in_row, *low_row, *high_row;

	for (i = 0; i < r; ++i) {
		in_row   = (float *)(in->data.ptr + i * in->step);
		low_row  = (float *)(low->data.ptr + i * low->step);
		high_row = (float *)(high->data.ptr + i * high->step);

		for (j = 0; j < c / 2; ++j) {
			a           = in_row[2 * j];
			b           = in_row[2 * j + 1];
			low_row[j]  = (a + b) * inv_sqrt2;
			high_row[j] = (a - b) * inv_sqrt2;
		}
	}
}

static void
dwt_col_pass(const CvMat *low, const CvMat *high,
             CvMat *LL, CvMat *LH, CvMat *HL, CvMat *HH)
{
	const float inv_sqrt2 = 1.0f / sqrtf(2.0f);
	int         r         = LL->rows;
	int         c         = LL->cols;
	int         i, j;
	float       a, b, c_val, d;
	float      *low0, *low1, *high0, *high1;
	float      *LL_row, *LH_row, *HL_row, *HH_row;

	for (i = 0; i < r; ++i) {
		low0   = (float *)(low->data.ptr + (2 * i) * low->step);
		low1   = (float *)(low->data.ptr + (2 * i + 1) * low->step);
		high0  = (float *)(high->data.ptr + (2 * i) * high->step);
		high1  = (float *)(high->data.ptr + (2 * i + 1) * high->step);

		LL_row = (float *)(LL->data.ptr + i * LL->step);
		LH_row = (float *)(LH->data.ptr + i * LH->step);
		HL_row = (float *)(HL->data.ptr + i * HL->step);
		HH_row = (float *)(HH->data.ptr + i * HH->step);

		for (j = 0; j < c; ++j) {
			a         = low0[j];
			b         = low1[j];
			c_val     = high0[j];
			d         = high1[j];

			LL_row[j] = (a + b) * inv_sqrt2;
			HL_row[j] = (a - b) * inv_sqrt2;
			LH_row[j] = (c_val + d) * inv_sqrt2;
			HH_row[j] = (c_val - d) * inv_sqrt2;
		}
	}
}

static void
haar_dwt2d(struct ProcCtx *ctx, const CvMat *input)
{
	const CvMat *src;
	int          lvl;

	src = input;
	for (lvl = 0; lvl < ctx->levels; ++lvl) {
		dwt_row_pass(src, ctx->dwt[lvl].low, ctx->dwt[lvl].high);
		dwt_col_pass(ctx->dwt[lvl].low, ctx->dwt[lvl].high,
		             ctx->dwt[lvl].LL, ctx->dwt[lvl].LH,
		             ctx->dwt[lvl].HL, ctx->dwt[lvl].HH);
		src = ctx->dwt[lvl].LL;
	}
}

static void
adaptive_threshold_detail(const CvMat *detail, struct ThreshScratch *s, IplImage *binary,
                          float base_thresh_factor, int block_size, float C)
{
	double min_val, max_val;
	double scale, shift;
	double otsu;

	cvAbsDiffS(detail, s->abs_detail, cvScalarAll(0));

	min_val = 0.0;
	max_val = 0.0;
	cvMinMaxLoc(s->abs_detail, &min_val, &max_val, NULL, NULL, NULL);

	if (max_val > min_val) {
		scale = 255.0 / (max_val - min_val);
		shift = -255.0 * min_val / (max_val - min_val);
		cvConvertScale(s->abs_detail, s->normalized, scale, shift);
	} else {
		cvZero(s->normalized);
	}

	if (block_size % 2 == 0)
		block_size += 1;

	cvAdaptiveThreshold(s->normalized, binary, 255, CV_ADAPTIVE_THRESH_GAUSSIAN_C,
	                    CV_THRESH_BINARY, block_size, C);

	otsu = cvThreshold(s->normalized, s->global_binary, 0, 255,
	                   CV_THRESH_BINARY | CV_THRESH_OTSU);
	cvThreshold(s->normalized, s->global_binary, otsu * base_thresh_factor, 255,
	            CV_THRESH_BINARY);
	cvAnd(binary, s->global_binary, binary, NULL);
}

static void
smart_morphology(const IplImage *binary, struct MorphSet *m, IplImage *out)
{
	cvCopy(binary, out, NULL);
	cvMorphologyEx(out, out, m->temp, m->k_text, CV_MOP_CLOSE, 1);
	cvMorphologyEx(out, out, m->temp, m->k_clean, CV_MOP_OPEN, 1);
	cvDilate(out, out, m->k_text, 1);
}

static int
rect_intensity_score(const IplImage *luma, CvRect r, float *bright_ratio, float *contrast)
{
	const unsigned char *row;
	int                  x, y, v;
	int                  bright, total;
	int                  max_v, min_v;

	if (r.x < 0 || r.y < 0 || r.x + r.width > luma->width || r.y + r.height > luma->height)
		return 0;

	bright = 0;
	total  = 0;
	max_v  = 0;
	min_v  = 255;

	for (y = r.y; y < r.y + r.height; y += 2) {
		row = (const unsigned char *)(luma->imageData + y * luma->widthStep);
		for (x = r.x; x < r.x + r.width; x += 2) {
			v = row[x];
			if (v > max_v)
				max_v = v;
			if (v < min_v)
				min_v = v;
			if (v > BRIGHT_THRESHOLD)
				bright += 1;
			total += 1;
		}
	}

	if (total == 0)
		return 0;

	*bright_ratio = (float)bright / (float)total;
	*contrast     = (float)(max_v - min_v) / 255.0f;
	return 1;
}

static int
imax(int a, int b)
{
	return a > b ? a : b;
}

static int
region_passes_geometry(CvRect scaled, float aspect, float fullness,
                       int orig_w, int orig_h)
{
	float center_x       = scaled.x + scaled.width / 2.0f;
	float frame_center_x = orig_w / 2.0f;
	int   bottom_limit   = (int)(orig_h * REGION_BOTTOM_LIMIT_PCT);
	int   top_limit      = (int)(orig_h * REGION_TOP_LIMIT_PCT);
	int   y_bottom       = scaled.y + scaled.height;
	int   min_w          = imax(REGION_MIN_WIDTH_FLOOR, (int)(orig_w * REGION_MIN_WIDTH_PCT));
	int   max_w          = (int)(orig_w * REGION_MAX_WIDTH_PCT);
	int   min_h          = imax(REGION_MIN_HEIGHT_FLOOR, (int)(orig_h * REGION_MIN_HEIGHT_PCT));
	int   max_h          = imax(REGION_MAX_HEIGHT_FLOOR, (int)(orig_h * REGION_MAX_HEIGHT_PCT));

	if (fabsf(center_x - frame_center_x) > frame_center_x * REGION_CENTER_TOL_PCT)
		return 0;

	return scaled.width  > min_w && scaled.width  < max_w
	    && scaled.height > min_h && scaled.height < max_h
	    && aspect        > REGION_MIN_ASPECT && aspect < REGION_MAX_ASPECT
	    && scaled.y      > top_limit
	    && y_bottom      < bottom_limit
	    && fullness      > REGION_MIN_FULLNESS;
}

static void
merge_horizontal_lines(struct TextRegionList *list, int max_gap_px)
{
	int    merged_any;
	int    i, j;
	int    yi_c, yj_c, hi, hj;
	int    gap_lr, gap_rl;
	int    x1, y1, x2, y2;
	int    irx2, iry2, jrx2, jry2;
	CvRect ri, rj;

	do {
		merged_any = 0;
		for (i = 0; i < list->count; ++i) {
			ri = list->regions[i].rect;
			for (j = i + 1; j < list->count; ++j) {
				rj   = list->regions[j].rect;
				yi_c = ri.y + ri.height / 2;
				yj_c = rj.y + rj.height / 2;
				hi   = ri.height;
				hj   = rj.height;

				if (abs(yi_c - yj_c) > (hi + hj) / 4)
					continue;
				if (hi > hj * 2 || hj > hi * 2)
					continue;

				gap_lr = rj.x - (ri.x + ri.width);
				gap_rl = ri.x - (rj.x + rj.width);
				if (gap_lr > max_gap_px && gap_rl > max_gap_px)
					continue;

				irx2                         = ri.x + ri.width;
				iry2                         = ri.y + ri.height;
				jrx2                         = rj.x + rj.width;
				jry2                         = rj.y + rj.height;

				x1                           = ri.x < rj.x ? ri.x : rj.x;
				y1                           = ri.y < rj.y ? ri.y : rj.y;
				x2                           = irx2 > jrx2 ? irx2 : jrx2;
				y2                           = iry2 > jry2 ? iry2 : jry2;

				list->regions[i].rect.x      = x1;
				list->regions[i].rect.y      = y1;
				list->regions[i].rect.width  = x2 - x1;
				list->regions[i].rect.height = y2 - y1;
				if (list->regions[j].confidence > list->regions[i].confidence)
					list->regions[i].confidence = list->regions[j].confidence;

				list->regions[j]  = list->regions[list->count - 1];
				list->count      -= 1;
				merged_any        = 1;
				j                 = i;
				ri                = list->regions[i].rect;
			}
		}
	} while (merged_any);
}

static void
nms(struct TextRegionList *list, float iou_thresh)
{
	int                   suppressed[MAX_REGIONS];
	int                   i, j, loser;
	struct TextRegionList out;

	memset(suppressed, 0, sizeof(suppressed));
	for (i = 0; i < list->count; ++i) {
		if (suppressed[i])
			continue;
		for (j = i + 1; j < list->count; ++j) {
			if (suppressed[j])
				continue;
			if (rect_iou(list->regions[i].rect, list->regions[j].rect) > iou_thresh) {
				loser             = list->regions[i].confidence
				                                    >= list->regions[j].confidence ?
				                            j :
				                            i;
				suppressed[loser] = 1;
				if (loser == i)
					break;
			}
		}
	}

	memset(&out, 0, sizeof(out));
	for (i = 0; i < list->count; ++i) {
		if (suppressed[i])
			continue;
		if (out.count >= MAX_REGIONS)
			break;
		out.regions[out.count++] = list->regions[i];
	}
	*list = out;
}

static void
detect_text_regions(struct ProcCtx *ctx, const IplImage *edge_map, const IplImage *luma,
                    int orig_h, int orig_w, int scale_factor,
                    struct TextRegionList *output)
{
	CvSeq                *contours;
	CvSeq                *c;
	CvRect                rect, scaled;
	struct TextRegionList raw;
	float                 aspect, contour_area, rect_area, fullness;
	float                 bright_ratio, contrast;
	float                 confidence;

	int                   bottom_limit = (int)(orig_h * REGION_BOTTOM_LIMIT_PCT);
	int                   top_limit    = (int)(orig_h * REGION_TOP_LIMIT_PCT);
	int                   bottom_band  = (int)(orig_h * REGION_BOTTOM_BAND_PCT);
	int                   min_h        = imax(REGION_MIN_HEIGHT_FLOOR,
	                                          (int)(orig_h * REGION_MIN_HEIGHT_PCT));
	int                   max_h        = imax(REGION_MAX_HEIGHT_FLOOR,
	                                          (int)(orig_h * REGION_MAX_HEIGHT_PCT));
	int                   y_bottom;

	cvCopy(edge_map, ctx->edge_copy, NULL);
	cvClearMemStorage(ctx->storage);
	contours = NULL;
	cvFindContours(ctx->edge_copy, ctx->storage, &contours, sizeof(CvContour),
	               CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));

	memset(&raw, 0, sizeof(raw));

	for (c = contours; c != NULL; c = c->h_next) {
		if (c->total < REGION_MIN_CONTOUR_LEN)
			continue;

		rect = cvBoundingRect(c, 0);
		if (rect.width <= 0 || rect.height <= 0)
			continue;

		scaled.x      = rect.x * scale_factor;
		scaled.y      = rect.y * scale_factor;
		scaled.width  = rect.width * scale_factor;
		scaled.height = rect.height * scale_factor;

		y_bottom      = scaled.y + scaled.height;
		if (scaled.y < top_limit || y_bottom > bottom_limit)
			continue;
		if (y_bottom < bottom_band)
			continue;
		if (scaled.height < min_h || scaled.height > max_h)
			continue;

		contour_area = (float)cvContourArea(c, CV_WHOLE_SEQ, 0);
		rect_area    = (float)rect.width * (float)rect.height;
		fullness     = rect_area > 0.0f ? contour_area / rect_area : 0.0f;
		if (fullness < REGION_MIN_FULLNESS)
			continue;

		if (raw.count < MAX_REGIONS) {
			raw.regions[raw.count].rect        = scaled;
			raw.regions[raw.count].confidence  = fullness;
			raw.regions[raw.count].is_stable   = 0;
			raw.count                         += 1;
		}
	}

	merge_horizontal_lines(&raw, imax(80, orig_w / 16));

	{
		struct TextRegionList filtered = { 0 };
		int                   k;
		for (k = 0; k < raw.count; ++k) {
			scaled = raw.regions[k].rect;
			aspect = (float)scaled.width / (float)scaled.height;
			if (!region_passes_geometry(scaled, aspect, raw.regions[k].confidence,
			                            orig_w, orig_h))
				continue;

			bright_ratio = 0.0f;
			contrast     = 0.0f;
			if (!rect_intensity_score(luma, scaled, &bright_ratio, &contrast))
				continue;
			if (bright_ratio < BRIGHT_RATIO_MIN || bright_ratio > BRIGHT_RATIO_MAX)
				continue;
			if (contrast < CONTRAST_MIN)
				continue;

			confidence  = 1.0f;
			confidence *= fminf(aspect / 6.0f, 1.0f);
			confidence *= fminf(scaled.height / 32.0f, 1.0f);
			confidence *= raw.regions[k].confidence;
			confidence *= fminf(contrast * 1.5f, 1.0f);
			confidence *= fminf(bright_ratio * 6.0f, 1.0f);

			if (filtered.count < MAX_REGIONS) {
				filtered.regions[filtered.count].rect        = scaled;
				filtered.regions[filtered.count].confidence  = confidence;
				filtered.regions[filtered.count].is_stable   = 0;
				filtered.count                              += 1;
			}
		}

		nms(&filtered, NMS_IOU_THRESH);
		*output = filtered;
	}
}

static int
tracker_alloc_slot(struct Tracker *tk)
{
	int t;

	for (t = 0; t < MAX_TRACKS; ++t) {
		if (t >= tk->count) {
			tk->count = t + 1;
			return t;
		}
		if (!tk->tracks[t].active)
			return t;
	}
	return -1;
}

static int
tracker_update(struct Tracker *tk, const struct TextRegionList *detections,
               struct TextRegionList *stable_out)
{
	int                   matched_det[MAX_REGIONS];
	int                   t, d, slot;
	int                   best;
	int                   new_captures;
	float                 best_iou, iou;
	float                 keep, kc;
	CvRect                dr, tr;
	struct TextRegionList out;

	memset(matched_det, 0, sizeof(matched_det));

	new_captures = 0;
	keep         = TRACK_EMA_KEEP;
	kc           = TRACK_CONF_EMA_KEEP;

	for (t = 0; t < tk->count; ++t) {
		if (!tk->tracks[t].active)
			continue;

		best     = -1;
		best_iou = TRACK_IOU_MIN;
		for (d = 0; d < detections->count; ++d) {
			if (matched_det[d])
				continue;
			iou = rect_iou(tk->tracks[t].rect, detections->regions[d].rect);
			if (iou > best_iou) {
				best_iou = iou;
				best     = d;
			}
		}

		if (best >= 0) {
			dr                         = detections->regions[best].rect;
			tr                         = tk->tracks[t].rect;
			tk->tracks[t].rect.x       = (int)(keep * tr.x + (1 - keep) * dr.x);
			tk->tracks[t].rect.y       = (int)(keep * tr.y + (1 - keep) * dr.y);
			tk->tracks[t].rect.width   = (int)(keep * tr.width + (1 - keep) * dr.width);
			tk->tracks[t].rect.height  = (int)(keep * tr.height + (1 - keep) * dr.height);
			tk->tracks[t].confidence   = kc * tk->tracks[t].confidence
			                             + (1 - kc) * detections->regions[best].confidence;
			tk->tracks[t].hits        += 1;
			tk->tracks[t].misses       = 0;
			matched_det[best]          = 1;
		} else {
			tk->tracks[t].misses += 1;
			if (tk->tracks[t].misses > TRACK_MAX_MISSES)
				tk->tracks[t].active = 0;
		}
	}

	for (d = 0; d < detections->count; ++d) {
		if (matched_det[d])
			continue;
		slot = tracker_alloc_slot(tk);
		if (slot < 0)
			continue;

		tk->tracks[slot].rect       = detections->regions[d].rect;
		tk->tracks[slot].confidence = detections->regions[d].confidence;
		tk->tracks[slot].hits       = 1;
		tk->tracks[slot].misses     = 0;
		tk->tracks[slot].active     = 1;
		tk->tracks[slot].captured   = 0;
	}

	memset(&out, 0, sizeof(out));
	for (t = 0; t < tk->count; ++t) {
		if (!tk->tracks[t].active)
			continue;
		if (tk->tracks[t].hits < TRACK_HIT_EMIT)
			continue;
		if (out.count >= MAX_REGIONS)
			break;

		out.regions[out.count].rect        = tk->tracks[t].rect;
		out.regions[out.count].confidence  = tk->tracks[t].confidence;
		out.regions[out.count].is_stable   = tk->tracks[t].hits >= TRACK_HIT_STABLE;
		out.count                         += 1;

		if (!tk->tracks[t].captured && tk->tracks[t].hits >= TRACK_HIT_STABLE) {
			tk->tracks[t].captured  = 1;
			new_captures           += 1;
		}
	}
	*stable_out = out;
	return new_captures;
}

static void
enhance_subtitles(const IplImage *frame, const struct TextRegionList *regions,
                  IplImage *gray_scratch, IplImage *bin_scratch, IplImage *output)
{
	CvRect         rect;
	int            idx, x, y;
	unsigned char *mask_row;
	unsigned char *out_row;
	unsigned char  val;

	cvCopy(frame, output, NULL);

	for (idx = 0; idx < regions->count; ++idx) {
		rect = regions->regions[idx].rect;
		if (rect.x < 0)
			rect.x = 0;
		if (rect.y < 0)
			rect.y = 0;
		if (rect.x + rect.width > frame->width)
			rect.width = frame->width - rect.x;
		if (rect.y + rect.height > frame->height)
			rect.height = frame->height - rect.y;

		cvSetImageROI((IplImage *)frame, rect);
		cvSetImageROI(gray_scratch, rect);
		cvCvtColor(frame, gray_scratch, CV_BGR2GRAY);
		cvEqualizeHist(gray_scratch, gray_scratch);

		cvSetImageROI(bin_scratch, rect);
		cvAdaptiveThreshold(gray_scratch, bin_scratch, 255, CV_ADAPTIVE_THRESH_GAUSSIAN_C,
		                    CV_THRESH_BINARY, 11, 2);

		for (y = 0; y < rect.height; ++y) {
			mask_row = (unsigned char *)(bin_scratch->imageData
			                             + (rect.y + y) * bin_scratch->widthStep
			                             + rect.x);
			out_row  = (unsigned char *)(output->imageData
			                             + (rect.y + y) * output->widthStep
			                             + rect.x * output->nChannels);
			for (x = 0; x < rect.width; ++x) {
				val                = mask_row[x] > 128 ? 255 : 0;
				out_row[x * 3]     = val;
				out_row[x * 3 + 1] = val;
				out_row[x * 3 + 2] = val;
			}
		}

		cvResetImageROI((IplImage *)frame);
		cvResetImageROI(gray_scratch);
		cvResetImageROI(bin_scratch);
	}
}

static int
procctx_init(struct ProcCtx *ctx, int W, int H, int levels)
{
	int    lvl;
	int    r, c;
	CvSize detail_size;

	memset(ctx, 0, sizeof(*ctx));
	ctx->W      = W;
	ctx->H      = H;
	ctx->levels = levels;
	ctx->dW     = W >> levels;
	ctx->dH     = H >> levels;
	detail_size = cvSize(ctx->dW, ctx->dH);

	ctx->gray   = cvCreateImage(cvSize(W, H), IPL_DEPTH_8U, 1);
	ctx->gray32 = cvCreateMat(H, W, CV_32F);

	ctx->dwt    = (struct DwtLevel *)calloc(levels, sizeof(struct DwtLevel));
	if (!ctx->dwt)
		return 0;

	for (lvl = 0; lvl < levels; ++lvl) {
		r                  = H >> lvl;
		c                  = W >> lvl;
		ctx->dwt[lvl].low  = cvCreateMat(r, c / 2, CV_32F);
		ctx->dwt[lvl].high = cvCreateMat(r, c / 2, CV_32F);
		ctx->dwt[lvl].LL   = cvCreateMat(r / 2, c / 2, CV_32F);
		ctx->dwt[lvl].LH   = cvCreateMat(r / 2, c / 2, CV_32F);
		ctx->dwt[lvl].HL   = cvCreateMat(r / 2, c / 2, CV_32F);
		ctx->dwt[lvl].HH   = cvCreateMat(r / 2, c / 2, CV_32F);
	}

	ctx->t_scratch.abs_detail    = cvCreateMat(ctx->dH, ctx->dW, CV_32F);
	ctx->t_scratch.normalized    = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->t_scratch.global_binary = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);

	ctx->binLH                   = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->binHL                   = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->binHH                   = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);

	ctx->morphLH                 = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->morphHL                 = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->morphHH                 = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);

	ctx->morphLH_set.temp        = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->morphLH_set.k_text      = cvCreateStructuringElementEx(
	        MORPH_LH_TEXT_W, MORPH_LH_TEXT_H,
	        MORPH_LH_TEXT_W / 2, MORPH_LH_TEXT_H / 2, CV_SHAPE_RECT, NULL);
	ctx->morphLH_set.k_clean = cvCreateStructuringElementEx(
	        MORPH_LH_CLEAN_W, MORPH_LH_CLEAN_H,
	        MORPH_LH_CLEAN_W / 2, MORPH_LH_CLEAN_H / 2, CV_SHAPE_RECT, NULL);

	ctx->morphHL_set.temp   = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->morphHL_set.k_text = cvCreateStructuringElementEx(
	        MORPH_HL_TEXT_W, MORPH_HL_TEXT_H,
	        MORPH_HL_TEXT_W / 2, MORPH_HL_TEXT_H / 2, CV_SHAPE_RECT, NULL);
	ctx->morphHL_set.k_clean = cvCreateStructuringElementEx(
	        MORPH_HL_CLEAN_W, MORPH_HL_CLEAN_H,
	        MORPH_HL_CLEAN_W / 2, MORPH_HL_CLEAN_H / 2, CV_SHAPE_RECT, NULL);

	ctx->morphHH_set.temp   = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->morphHH_set.k_text = cvCreateStructuringElementEx(
	        MORPH_HH_TEXT_W, MORPH_HH_TEXT_H,
	        MORPH_HH_TEXT_W / 2, MORPH_HH_TEXT_H / 2, CV_SHAPE_RECT, NULL);
	ctx->morphHH_set.k_clean = cvCreateStructuringElementEx(
	        MORPH_HH_CLEAN_W, MORPH_HH_CLEAN_H,
	        MORPH_HH_CLEAN_W / 2, MORPH_HH_CLEAN_H / 2, CV_SHAPE_RECT, NULL);

	ctx->edge_map     = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->edge_copy    = cvCreateImage(detail_size, IPL_DEPTH_8U, 1);
	ctx->line_kernel  = cvCreateStructuringElementEx(LINE_KERNEL_W, LINE_KERNEL_H,
	                                                 LINE_KERNEL_W / 2, LINE_KERNEL_H / 2,
	                                                 CV_SHAPE_RECT, NULL);
	ctx->clean_kernel = cvCreateStructuringElementEx(CLEAN_KERNEL_W, CLEAN_KERNEL_H,
	                                                 CLEAN_KERNEL_W / 2, CLEAN_KERNEL_H / 2,
	                                                 CV_SHAPE_RECT, NULL);
	ctx->storage      = cvCreateMemStorage(0);

	ctx->display      = cvCreateImage(cvSize(W, H), IPL_DEPTH_8U, 3);
	ctx->enh_gray     = cvCreateImage(cvSize(W, H), IPL_DEPTH_8U, 1);
	ctx->enh_binary   = cvCreateImage(cvSize(W, H), IPL_DEPTH_8U, 1);

	ctx->vchan        = cvCreateImage(cvSize(W, H), IPL_DEPTH_8U, 1);
	ctx->salience     = cvCreateImage(cvSize(W, H), IPL_DEPTH_8U, 1);

	ctx->ppm_buf_sz   = (size_t)W * (size_t)H * 3;
	ctx->ppm_buf      = (unsigned char *)malloc(ctx->ppm_buf_sz);

	return 1;
}

static void
procctx_free(struct ProcCtx *ctx)
{
	int lvl;

	safe_release_image(&ctx->gray);
	safe_release_mat(&ctx->gray32);

	if (ctx->dwt) {
		for (lvl = 0; lvl < ctx->levels; ++lvl) {
			safe_release_mat(&ctx->dwt[lvl].low);
			safe_release_mat(&ctx->dwt[lvl].high);
			safe_release_mat(&ctx->dwt[lvl].LL);
			safe_release_mat(&ctx->dwt[lvl].LH);
			safe_release_mat(&ctx->dwt[lvl].HL);
			safe_release_mat(&ctx->dwt[lvl].HH);
		}
		free(ctx->dwt);
		ctx->dwt = NULL;
	}

	safe_release_mat(&ctx->t_scratch.abs_detail);
	safe_release_image(&ctx->t_scratch.normalized);
	safe_release_image(&ctx->t_scratch.global_binary);

	safe_release_image(&ctx->binLH);
	safe_release_image(&ctx->binHL);
	safe_release_image(&ctx->binHH);

	safe_release_image(&ctx->morphLH_set.temp);
	safe_release_kernel(&ctx->morphLH_set.k_text);
	safe_release_kernel(&ctx->morphLH_set.k_clean);
	safe_release_image(&ctx->morphHL_set.temp);
	safe_release_kernel(&ctx->morphHL_set.k_text);
	safe_release_kernel(&ctx->morphHL_set.k_clean);
	safe_release_image(&ctx->morphHH_set.temp);
	safe_release_kernel(&ctx->morphHH_set.k_text);
	safe_release_kernel(&ctx->morphHH_set.k_clean);

	safe_release_image(&ctx->morphLH);
	safe_release_image(&ctx->morphHL);
	safe_release_image(&ctx->morphHH);

	safe_release_image(&ctx->edge_map);
	safe_release_image(&ctx->edge_copy);
	safe_release_kernel(&ctx->line_kernel);
	safe_release_kernel(&ctx->clean_kernel);
	if (ctx->storage)
		cvReleaseMemStorage(&ctx->storage);

	safe_release_image(&ctx->display);
	safe_release_image(&ctx->enh_gray);
	safe_release_image(&ctx->enh_binary);

	safe_release_image(&ctx->vchan);
	safe_release_image(&ctx->salience);

	if (ctx->ppm_buf) {
		free(ctx->ppm_buf);
		ctx->ppm_buf = NULL;
	}
}

static void
compute_salience_and_v(const IplImage *frame, IplImage *salience, IplImage *vchan)
{
	int            x, y;
	int            W        = frame->width;
	int            H        = frame->height;
	int            src_step = frame->widthStep;
	int            sal_step = salience->widthStep;
	int            v_step   = vchan->widthStep;
	const uint8_t *src_row;
	uint8_t       *sal_row;
	uint8_t       *v_row;
	int            b, g, r, vmax, vmin, chroma, sal;

	for (y = 0; y < H; ++y) {
		src_row = (const uint8_t *)(frame->imageData + y * src_step);
		sal_row = (uint8_t *)(salience->imageData + y * sal_step);
		v_row   = (uint8_t *)(vchan->imageData + y * v_step);

		for (x = 0; x < W; ++x) {
			b    = src_row[x * 3];
			g    = src_row[x * 3 + 1];
			r    = src_row[x * 3 + 2];

			vmax = b > g ? b : g;
			if (r > vmax)
				vmax = r;
			vmin = b < g ? b : g;
			if (r < vmin)
				vmin = r;

			chroma = vmax - vmin;
			sal    = (vmax >> 1) + chroma;
			if (sal > 255)
				sal = 255;

			v_row[x]   = (uint8_t)vmax;
			sal_row[x] = (uint8_t)sal;
		}
	}
}

static int
process_frame(struct ProcCtx *ctx, const IplImage *frame, struct Tracker *tracker,
              float dwt_thresh, struct TextRegionList *regions_out)
{
	struct TextRegionList detections;
	int                   scale_factor;
	CvMat                *LH, *HL, *HH;

	compute_salience_and_v(frame, ctx->salience, ctx->vchan);

	cvCopy(ctx->salience, ctx->gray, NULL);
	cvEqualizeHist(ctx->gray, ctx->gray);
	cvSmooth(ctx->gray, ctx->gray, CV_MEDIAN, 3, 0, 0, 0);
	cvConvertScale(ctx->gray, ctx->gray32, 1.0 / 255.0, 0);

	haar_dwt2d(ctx, ctx->gray32);
	LH = ctx->dwt[ctx->levels - 1].LH;
	HL = ctx->dwt[ctx->levels - 1].HL;
	HH = ctx->dwt[ctx->levels - 1].HH;

	adaptive_threshold_detail(LH, &ctx->t_scratch, ctx->binLH,
	                          dwt_thresh, THRESH_BLOCK_SIZE, THRESH_C);
	adaptive_threshold_detail(HL, &ctx->t_scratch, ctx->binHL,
	                          dwt_thresh, THRESH_BLOCK_SIZE, THRESH_C);
	adaptive_threshold_detail(HH, &ctx->t_scratch, ctx->binHH,
	                          dwt_thresh * THRESH_HH_FACTOR, THRESH_BLOCK_SIZE, THRESH_C);

	smart_morphology(ctx->binLH, &ctx->morphLH_set, ctx->morphLH);
	smart_morphology(ctx->binHL, &ctx->morphHL_set, ctx->morphHL);
	smart_morphology(ctx->binHH, &ctx->morphHH_set, ctx->morphHH);

	cvAnd(ctx->morphLH, ctx->morphHL, ctx->edge_map, NULL);
	cvAnd(ctx->edge_map, ctx->morphHH, ctx->edge_map, NULL);

	cvDilate(ctx->edge_map, ctx->edge_map, ctx->line_kernel, LINE_DILATE_ITERS);
	cvMorphologyEx(ctx->edge_map, ctx->edge_map, NULL, ctx->line_kernel, CV_MOP_CLOSE, 1);
	cvErode(ctx->edge_map, ctx->edge_map, ctx->clean_kernel, 1);

	scale_factor = 1 << ctx->levels;
	memset(&detections, 0, sizeof(detections));
	detect_text_regions(ctx, ctx->edge_map, ctx->vchan, frame->height, frame->width,
	                    scale_factor, &detections);

	return tracker_update(tracker, &detections, regions_out);
}

static int
save_ppm(const char *filename, const IplImage *img, unsigned char *swap_buf)
{
	FILE                *f;
	const unsigned char *src;
	unsigned char       *dst;
	int                  x, y;
	int                  W, H;
	size_t               row_bytes;

	if (!img || img->nChannels < 3 || img->depth != IPL_DEPTH_8U || !swap_buf)
		return 0;

	W         = img->width;
	H         = img->height;
	row_bytes = (size_t)W * 3;

	for (y = 0; y < H; ++y) {
		src = (const unsigned char *)(img->imageData + y * img->widthStep);
		dst = swap_buf + (size_t)y * row_bytes;
		for (x = 0; x < W; ++x) {
			dst[x * 3]     = src[x * 3 + 2];
			dst[x * 3 + 1] = src[x * 3 + 1];
			dst[x * 3 + 2] = src[x * 3];
		}
	}

	f = fopen(filename, "wb");
	if (!f)
		return 0;

	fprintf(f, "P6\n%d %d\n255\n", W, H);
	fwrite(swap_buf, 1, row_bytes * (size_t)H, f);
	fclose(f);
	return 1;
}

static void
draw_overlay(IplImage *display, const struct TextRegionList *regions,
             float min_confidence, const CvFont *font,
             int frame_count, double pos_s, double total_s, double fps)
{
	struct TextRegion r;
	char              conf_text[32];
	char              info[160];
	int               i, green, red;

	for (i = 0; i < regions->count; ++i) {
		r = regions->regions[i];
		if (r.confidence < min_confidence)
			continue;

		green = (int)(r.confidence * 255.0f);
		if (green > 255)
			green = 255;
		red = 255 - green;

		cvRectangle(display,
		            cvPoint(r.rect.x, r.rect.y),
		            cvPoint(r.rect.x + r.rect.width, r.rect.y + r.rect.height),
		            CV_RGB(0, green, red), r.is_stable ? 3 : 1, 8, 0);

		snprintf(conf_text, sizeof(conf_text), "C: %.2f", r.confidence);
		cvPutText(display, conf_text, cvPoint(r.rect.x, r.rect.y - 5),
		          font, CV_RGB(0, 255, 0));
	}

	snprintf(info, sizeof(info),
	         "Frame: %d | Regions: %d | FPS: %.1f | %.1f/%.1fs",
	         frame_count, regions->count, fps, pos_s, total_s);
	cvPutText(display, info, cvPoint(10, 30), font, CV_RGB(0, 255, 255));
}

int
main(int argc, char **argv)
{
	struct Video          v;
	struct ProcCtx        ctx;
	struct Tracker        tracker;
	struct TextRegionList regions;
	CvFont                font;
	IplImage             *frame;
	const char           *win_name = "Subtitle Extraction";
	int                   dwt_thresh_int;
	int                   min_conf_int;
	int                   pos_int;
	int                   shown_pos;
	int                   paused;
	int                   enable_color_enhancement;
	int                   frame_count;
	int                   processed_frames;
	int                   new_caps;
	int                   key;
	float                 dwt_threshold_factor;
	float                 min_confidence;
	double                processing_sum;
	double                frame_start, frame_end;
	double                avg_time;
	char                  filename[128];

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <video_path>\n", argv[0]);
		return 1;
	}

	if (!video_open(&v, argv[1])) {
		fprintf(stderr, "Error opening video file: %s\n", argv[1]);
		video_close(&v);
		return 1;
	}

	if (!procctx_init(&ctx, v.width, v.height, DWT_LEVELS)) {
		fprintf(stderr, "Failed to allocate processing context\n");
		video_close(&v);
		return 1;
	}

	ensure_dir("screenshots");

	dwt_thresh_int           = INIT_DWT_THRESH_INT;
	min_conf_int             = INIT_MIN_CONF_INT;
	pos_int                  = 0;
	shown_pos                = 0;
	paused                   = 0;
	enable_color_enhancement = 1;

	cvNamedWindow(win_name, CV_WINDOW_NORMAL);
	cvResizeWindow(win_name, WIN_W, WIN_H);
	cvCreateTrackbar("DWT Thresh x100", win_name, &dwt_thresh_int, 100, NULL);
	cvCreateTrackbar("Min Confidence x100", win_name, &min_conf_int, 100, NULL);
	cvCreateTrackbar("Position", win_name, &pos_int, POSITION_TICKS, NULL);

	cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5, 0, 1, 8);

	frame            = cvCreateImage(cvSize(v.width, v.height), IPL_DEPTH_8U, 3);
	frame_count      = 0;
	processed_frames = 0;
	processing_sum   = 0.0;
	memset(&tracker, 0, sizeof(tracker));

	while (1) {
		if (pos_int != shown_pos && v.duration_us > 0) {
			int64_t target_us = (int64_t)pos_int * v.duration_us / POSITION_TICKS;
			seek_reset(&v, &tracker, target_us);
			shown_pos = pos_int;
		}

		if (paused) {
			key = cvWaitKey(30);
			goto handle_key;
		}

		if (!video_read(&v, frame))
			break;

		frame_start  = now_seconds();
		frame_count += 1;
		if (frame_count % FRAME_SKIP_MOD != 0)
			continue;

		dwt_threshold_factor = dwt_thresh_int / 100.0f;
		min_confidence       = min_conf_int / 100.0f;

		memset(&regions, 0, sizeof(regions));
		new_caps = process_frame(&ctx, frame, &tracker, dwt_threshold_factor, &regions);

		if (enable_color_enhancement && regions.count > 0) {
			enhance_subtitles(frame, &regions, ctx.enh_gray, ctx.enh_binary, ctx.display);
		} else {
			cvCopy(frame, ctx.display, NULL);
		}

		{
			double pos_s   = v.cur_pts_us / 1e6;
			double total_s = v.duration_us > 0 ? v.duration_us / 1e6 : 0.0;
			draw_overlay(ctx.display, &regions, min_confidence, &font,
			             frame_count, pos_s, total_s, v.fps);
		}

		if (new_caps > 0) {
			snprintf(filename, sizeof(filename),
			         "screenshots/frame_%d.ppm", frame_count);
			if (save_ppm(filename, ctx.display, ctx.ppm_buf))
				printf("Saved screenshot: %s (+%d new)\n", filename, new_caps);
		}

		cvShowImage(win_name, ctx.display);

		if (v.duration_us > 0) {
			int new_pos = (int)((v.cur_pts_us * POSITION_TICKS) / v.duration_us);
			if (new_pos < 0)
				new_pos = 0;
			if (new_pos > POSITION_TICKS)
				new_pos = POSITION_TICKS;
			if (new_pos != pos_int) {
				cvSetTrackbarPos("Position", win_name, new_pos);
				shown_pos = new_pos;
			} else {
				shown_pos = pos_int;
			}
		}

		frame_end         = now_seconds();
		processing_sum   += (frame_end - frame_start) * 1000.0;
		processed_frames += 1;
		if (processed_frames % 30 == 0) {
			avg_time = processing_sum / processed_frames;
			printf("Avg processing time: %.2f ms\n", avg_time);
		}

		key = cvWaitKey(1);

handle_key:
		if (key == 27)
			break;
		if (key == 'c')
			enable_color_enhancement = !enable_color_enhancement;
		if (key == ' ' || key == 'p')
			paused = !paused;
		if (key == 'a' || key == 'h')
			seek_reset(&v, &tracker, v.cur_pts_us - SEEK_STEP_US);
		if (key == 'd' || key == 'l')
			seek_reset(&v, &tracker, v.cur_pts_us + SEEK_STEP_US);
		if (key == 'A' || key == 'H')
			seek_reset(&v, &tracker, v.cur_pts_us - SEEK_BIG_FACTOR * SEEK_STEP_US);
		if (key == 'D' || key == 'L')
			seek_reset(&v, &tracker, v.cur_pts_us + SEEK_BIG_FACTOR * SEEK_STEP_US);
	}

	safe_release_image(&frame);
	procctx_free(&ctx);
	video_close(&v);
	cvDestroyAllWindows();

	return 0;
}
