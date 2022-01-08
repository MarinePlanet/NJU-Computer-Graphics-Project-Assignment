#!/usr/bin/env python
# -*- coding:utf-8 -*-

# 本文件只允许依赖math库
import math


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'，此处的'Naive'仅作为示例，测试时不会出现
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    # print(x0, y0, x1, y1)
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, int(y0 + k * (x - x0))))




    elif algorithm == 'DDA':
        if x0 == x1:
            if y0 > y1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            for y in range(y0, y1+1, 1):
                result.append((x0, y))
        elif y0 == y1:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            for x in range(x0, x1+1, 1):
                result.append((x, y0))
        else:
            k = float((y1-y0))/(x1-x0)
            if abs(k) > 1:
                if k > 0:
                    if y0 > y1:
                        x0, y0, x1, y1 = x1, y1, x0, y0
                    x = x0
                    result.append((x0, y0))
                    for y in range(y0+1, y1+1):
                        x += float(1/k)
                        result.append((round(x), y))
                else:
                    if y0 <= y1:
                        x0, y0, x1, y1 = x1, y1, x0, y0
                    x = x0
                    result.append((x0, y0))
                    for y in range(y0-1, y1-1, -1):
                        x -= float(1/k)
                        result.append((round(x), y))
            elif abs(k) < 1:
                if k > 0:
                    if x0 > x1:
                        x0, y0, x1, y1 = x1, y1, x0, y0
                    y = y0
                    result.append((x0, y0))
                    for x in range(x0+1, x1+1):
                        y += k
                        result.append((x, round(y)))
                else:
                    if x0 > x1:
                        x0, y0, x1, y1 = x1, y1, x0, y0
                    y = y0
                    result.append((x0, y0))
                    for x in range(x0+1, x1+1):
                        y += k
                        result.append((x, round(y)))
            else:
                if k == 1:
                    if x0 > x1:
                        x0, y0, x1, y1 = x1, y1, x0, y0
                    result.append((x0, y0))
                    y = y0
                    for x in range(x0+1, x1+1):
                        y += 1
                        result.append((x, y))
                else:
                    if x0 > x1:
                        x0, y0, x1, y1 = x1, y1, x0, y0
                    result.append((x0, y0))
                    y = y0
                    for x in range(x0+1, x1+1):
                        y -= 1
                        result.append((x, y))





    elif algorithm == 'Bresenham':
        if x0 == x1:
            if y0 > y1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            for y in range(y0, y1+1, 1):
                result.append((x0, y))
        elif y0 == y1:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            for x in range(x0, x1+1, 1):
                result.append((x, y0))
        else:
            k = float((y1-y0))/(x1-x0)
            # print(f'k={k}')
            if (k > 0) and (k < 1):
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                delta_x = x1 - x0
                delta_y = y1 - y0
                incre_1 = 2*delta_y
                incre_2 = 2*delta_y - 2*delta_x
                p = 2*delta_y - delta_x
                x = x0
                y = y0
                while x <= x1:
                    result.append((x, y))
                    if p >= 0:
                        x += 1
                        y += 1
                        p += incre_2
                    else:
                        x += 1
                        p += incre_1

            elif k > 1:
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                delta_x = x1 - x0
                delta_y = y1 - y0
                incre_1 = 2*delta_x
                incre_2 = 2*delta_x - 2*delta_y
                p = 2*delta_x - delta_y
                x = x0
                y = y0
                while y <= y1:
                    result.append((x, y))
                    if p < 0:
                        y += 1
                        p += incre_1
                    else:
                        x += 1
                        y += 1
                        p += incre_2

            elif (k < 0) and (k > -1):
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                delta_x = x1 - x0
                delta_y = y1 - y0
                incre_1 = (-2)*delta_y
                incre_2 = (-2)*delta_x - 2*delta_y
                p = (-2)*delta_y - delta_x
                x = x0
                y = y0
                while x<=x1:
                    result.append((x, y))
                    if p < 0:
                        x += 1
                        p += incre_1
                    else:
                        x += 1
                        y -= 1
                        p += incre_2

            elif k < -1:
                if y0 > y1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                delta_x = x1 - x0
                delta_y = y1 - y0
                incre_1 = (-2)*delta_x
                incre_2 = (-2)*delta_x - 2*delta_y
                p = (-2)*delta_x - delta_y
                x = x0
                y = y0
                while y <= y1:
                    result.append((x, y))
                    if p < 0 :
                        y += 1
                        p += incre_1
                    else:
                        x -= 1
                        y += 1
                        p += incre_2

            elif k == 1:
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                y = y0
                result.append((x0, y))
                for x in range(x0+1, x1+1):
                    y += 1
                    result.append((x, y))

            else:
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                y = y0
                result.append((x0, y))
                for x in range(x0+1, x1+1):
                    y -= 1
                    result.append((x, y))
    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result


def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    xc = (float(x0) + float(x1))/2
    yc = (float(y0) + float(y1))/2
    rx = abs(float(x1) - xc)
    ry = abs(float(y0) - yc)
    p1 = (ry*ry)-(rx*rx*ry)+float(rx*rx/4)
    # print(f'xc={xc},yc={yc}')
    # print(f'rx={rx},ry={ry}')
    result.append((round(xc), round(yc + ry)))
    result.append((round(xc), round(yc - ry)))
    result.append((round(xc + rx), round(yc)))
    result.append((round(xc - rx), round(yc)))
    x = 0
    y = ry
    while (2*ry*ry*x) < (2*rx*rx*y):
        if p1 < 0:
            p1 += (2*ry*ry*x + 3*ry*ry)
            x += 1
            result.append((round(xc + x), round(yc + y)))
            result.append((round(xc + x), round(yc - y)))
            result.append((round(xc - x), round(yc - y)))
            result.append((round(xc - x), round(yc + y)))
        else:
            p1 += (2*ry*ry*x + 3*ry*ry - 2*rx*rx*y + 2*rx*rx)
            x += 1
            y -= 1
            result.append((round(xc + x), round(yc + y)))
            result.append((round(xc + x), round(yc - y)))
            result.append((round(xc - x), round(yc - y)))
            result.append((round(xc - x), round(yc + y)))
    # print(f'x={x},y={y}')
    # x2, y2 = result[-4]
    x2 = float(x)
    y2 = float(y)
    p2 = ry*ry*(x2+0.5)*(x2+0.5)+rx*rx*(y2-1)*(y2-1)-rx*rx*ry*ry
    # print(f'p2={p2}')
    # x = float(x2)
    # y = float(y2)
    # p2 = ry * ry * (x + 0.5) * (x + 0.5) + rx * rx * (y - 1) - rx * rx * ry * ry
    while y > 0:
        if p2 > 0:
            # print(f'incre_1:p2={p2}')
            p2 += (3*rx*rx - 2*rx*rx*y)
            y -= 1
            result.append((round(xc + x), round(yc + y)))
            result.append((round(xc + x), round(yc - y)))
            result.append((round(xc - x), round(yc - y)))
            result.append((round(xc - x), round(yc + y)))
            # print(f'{result[-4]}')
        else:
            # print(f'incre_2:p2={p2}')
            p2 += (2*ry*ry*x + 2*ry*ry + 3*rx*rx - 2*rx*rx*y)
            x += 1
            y -= 1
            result.append((round(xc + x), round(yc + y)))
            result.append((round(xc + x), round(yc - y)))
            result.append((round(xc - x), round(yc - y)))
            result.append((round(xc - x), round(yc + y)))
            # print(f'{result[-4]}')

    return result


def deCasteljau(n, i, u, p_list):
    if n == 0:
        x, y = p_list[i]
        return (x, y)
    else:
        x1, y1 = deCasteljau(n-1, i, u, p_list)
        x2, y2 = deCasteljau(n-1, i+1, u, p_list)
        x = (1-u)*x1 + u*x2
        y = (1-u)*y1 + u*y2
        return (x, y)


def de_boor_cox(u, k, i):
    if k == 1:
        if (u >= i) and (u < i + 1):
            return 1
        else:
            return 0
    elif k >= 2:
        c1 = (u - i)/(i + k - i - 1)
        c2 = (i + k - u)/(i + k - i - 1)
        return c1 * de_boor_cox(u, k-1, i)+c2 * de_boor_cox(u, k-1, i+1)


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    if len(p_list) == 1:
        return p_list
    if algorithm == 'Bezier':
        n = len(p_list) - 1
        for i in range(1, 10001):
            u = float(i/10000)
            x, y = deCasteljau(n, 0, u, p_list)
            result.append((round(x), round(y)))
    elif algorithm == 'B-spline':
        n = len(p_list)
        k = 4
        u = k - 1
        while u < n:
            x, y = 0, 0
            for i in range(0, n):
                dbc = de_boor_cox(u, k, i)
                xi, yi = p_list[i]
                x += xi*dbc
                y += yi*dbc
            result.append((round(x), round(y)))
            u += 1/1000

    return result


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    result = []
    for i in range(len(p_list)):
        xi, yi = p_list[i]
        result.append((xi+dx, yi+dy))
    return result


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    result = []
    for i in range(len(p_list)):
        xi, yi = p_list[i]
        dx = xi - x
        dy = yi - y
        theta = 2*math.pi*(float(r/360))
        new_dx = dx * math.cos(theta) - dy * math.sin(theta)
        new_dy = dx * math.sin(theta) + dy * math.cos(theta)
        result.append((round(x + new_dx), round(y + new_dy)))
    return result


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    result = []
    for i in range(len(p_list)):
        xi, yi = p_list[i]
        dx = xi - x
        dy = yi - y
        result.append((round(x + s * dx), round(y + s * dy)))
    return result


def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    # pass
    result = []
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    # print(f'x0,y0={p_list[0]}, x1,y1={p_list[1]}')
    if algorithm == 'Cohen-Sutherland':
        # print(f'Cohen-Sutherland p0 ={p_list[0]}, p1 ={p_list[1]}')
        # print(f'x_min={x_min},x_max={x_max},y_min={y_min},y_max={y_max}')
        while True:
            p0_byte = 0
            p1_byte = 0
            if x0 < x_min:
                p0_byte += 1
            if x0 > x_max:
                p0_byte += 2
            if y0 < y_min:
                p0_byte += 4
            if y0 > y_max:
                p0_byte += 8

            if x1 < x_min:
                p1_byte += 1
            if x1 > x_max:
                p1_byte += 2
            if y1 < y_min:
                p1_byte += 4
            if y1 > y_max:
                p1_byte += 8

            # print(f'p0,p1=({x0},{y0}),({x1},{y1})')
            # print(f'p0_byte={p0_byte},p1_byte={p1_byte}')

            if (p0_byte == 0) and (p1_byte == 0):
                result.append((x0, y0))
                result.append((x1, y1))
                break

            if p0_byte & p1_byte != 0:
                result.append((0, 0))
                result.append((0, 0))
                break

            else:
                if p0_byte == 0:
                    # print('Switched')
                    p0_byte, p1_byte = p1_byte, p0_byte
                    x0, y0, x1, y1 = x1, y1, x0, y0

                if p0_byte & 1 != 0:
                    y0 = round((x_min - x0) * (y1 - y0) / (x1 - x0) + y0)
                    x0 = x_min

                    # print(f'y1-y0={y1-y0},x1-x0={x1-x0}')
                    # print(f'test:{(y1-y0)/(x1-x0)}')
                    # print(f'case1:x0={x0},y0={y0}')

                if p0_byte & 2 != 0:
                    y0 = round((x_max - x0) * (float(y1) - y0) / (x1 - x0) + y0)
                    x0 = x_max


                if p0_byte & 4 != 0:
                    x0 = round((y_min-y0)*(float(x1)-x0)/(y1-y0) + x0)
                    y0 = y_min

                if p0_byte & 8 != 0:
                    x0 = round((y_max-y0)*(float(x1)-x0)/(y1-y0) + x0)
                    y0 = y_max
        # print(f'Cohen-Sutherland result = {result}')
        return result
    elif algorithm == 'Liang-Barsky':
        # result.append((0, 0))
        # result.append((0, 0))
        delta_x = x1 - x0
        delta_y = y1 - y0
        # print(f'delta_x = {delta_x}, delta_y={delta_y}')
        if delta_x == 0:
            if y0 > y1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            if (x0 >= x_min) and (x0 <= x_max) and (y0 <= y_max) and (y1 >= y_min):
                y_low = y0 if y0 > y_min else y_min
                y_high = y1 if y1 < y_max else y_max
                result.append((x0, y_low))
                result.append((x0, y_high))
            else:
                result.append((0, 0))
                result.append((0, 0))
        elif delta_y == 0:
            if x0 < x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            if (y0 >= y_min) and (y0 <= y_max) and (x0 < x_max) and (x0 > x_min):
                x_low = x0 if x0 > x_min else x_min
                x_high = x1 if x1 < x_max else x_max
                result.append((x_low, y0))
                result.append((x_high, y0))
            else:
                result.append((0, 0))
                result.append((0, 0))
        else:
            t_start = [0, 1, 2]
            t_end = [0, 1, 2]
            t_start[0] = 0.0
            t_end[0] = 1.0

            QL = -delta_x
            QR = delta_x
            QB = -delta_y
            QT = delta_y
            DL = x0 - x_min
            DR = x_max - x0
            DB = y0 - y_min
            DT = y_max - y0

            TL = float(DL)/QL
            TR = float(DR)/QR
            TB = float(DB)/QB
            TT = float(DT)/QT

            # print(f'TL={TL}')
            # print(f'TR={TR}')
            # print(f'TB={TB}')
            # print(f'TT={TT}')
            if delta_x >= 0:
                t_start[1] = TL
                t_end[1] = TR
            else:
                t_start[1] = TR
                t_end[1] = TL

            if delta_y >= 0:
                t_start[2] = TB
                t_end[2] = TT
            else:
                t_start[2] = TT
                t_end[2] = TB

            t0 = max(t_start)
            t1 = min(t_end)
            # print(f't0={t0},t1={t1}')
            if t0 < t1:
                result.append((round(x0 + t0*delta_x), round(y0 + t0*delta_y)))
                result.append((round(x0 + t1*delta_x), round(y0 + t1*delta_y)))
            # print(f' result={result}')
            else:
                result.append((0, 0))
                result.append((0, 0))
        return result
