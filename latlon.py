# -*- coding: utf-8 -*-
# Paul Reed
# 28 May 2019

import math

def sec(x): return 1/math.cos(x)

def getShift(x, y):
  e_idx = int(x/1000.0)
  n_idx = int(y/1000.0)
  x0 = e_idx * 1000
  y0 = n_idx * 1000
  shift_indexes = [
    [e_idx, n_idx], \
    [e_idx + 1, n_idx], \
    [e_idx + 1, n_idx + 1], \
    [e_idx, n_idx + 1]]
  r_nums = []
  s = []
  with open('OSTN15.csv', 'r') as file:
    lines = file.readlines()
    for i in range(4):
      r_num = shift_indexes[i][0] + (shift_indexes[i][1] * 701) + 1
      e_shift = float(lines[r_num].rstrip().split(",")[0])
      n_shift = float(lines[r_num].rstrip().split(",")[1])
      s.append([e_shift, n_shift])
    dx = x - x0
    dy = y - y0
    t = dx / 1000.0
    u = dy / 1000.0
    se = (1-t)*(1-u)*s[0][0]+t*(1-u)*s[1][0]+t*u*s[2][0]+(1-t)*u*s[3][0]
    sn = (1-t)*(1-u)*s[0][1]+t*(1-u)*s[1][1]+t*u*s[2][1]+(1-t)*u*s[3][1]
    return {'se': se, 'sn': sn}
  return None

def ETRS89ToOSGB36(latitude, longitude):
  lat = latitude * (math.pi/180)
  lon = longitude * (math.pi/180)
  a = 6378137.0000
  b = 6356752.3141
  e2 = (a**2-b**2)/a**2
  F0 = 0.9996012717
  lat0 = 49 * (math.pi/180)
  lon0 = -2 * (math.pi/180)
  E0 = 400000
  N0 = -100000

  n = (a-b)/(a+b)
  v = a*F0*(1-(e2*(math.sin(lat)**2)))**-0.5
  p = a*F0*(1-e2)*(1-(e2*(math.sin(lat)**2)))**-1.5
  n2 = v/p-1
  M = b*F0*((1+n+5.0/4.0*n**2+5.0/4.0*n**3)*(lat-lat0)-(3*n+3*n**2+21.0/8.0* \
    n**3)*math.sin(lat-lat0)*math.cos(lat+lat0)+(15.0/8.0*n**2+15.0/8.0*n**3 \
    ) * math.sin(2*(lat-lat0))*math.cos(2*(lat+lat0))-35.0/24.0*n**3* \
    math.sin(3*(lat-lat0))*math.cos(3*(lat+lat0)))
  I = M + N0
  II = v/2*math.sin(lat)*math.cos(lat)
  III = v/24*math.sin(lat)*math.cos(lat)**3*(5-math.tan(lat)**2+9*n2)
  IIIA = v/720*math.sin(lat)*math.cos(lat)**5*(61-58*math.tan(lat)**2+ \
    math.tan(lat)**4)
  IV = v * math.cos(lat)
  V = v/6*math.cos(lat)**3*(v/p-math.tan(lat)**2)
  VI = v/120*math.cos(lat)**5*(5-18*math.tan(lat)**2+math.tan(lat)**4+14*n2- \
    58*(math.tan(lat)**2*n2))

  N = I+II*(lon-lon0)**2+III*(lon-lon0)**4+IIIA*(lon-lon0)**6
  E = E0+IV*(lon-lon0)+V*(lon-lon0)**3+VI*(lon-lon0)**5

  shift = getShift(E, N)

  E = float('%.3f'%(E + shift['se']))
  N = float('%.3f'%(N + shift['sn']))

  return { 'EN':[E, N], 'LL':[latitude, longitude]}

def OSGB36ToETRS89(easting, northing):
  prev = getShift(easting, northing)
  new = getShift(easting - prev['se'], northing - prev['sn'])
  while abs(prev['se'] - new['se']) > 0.0001 or \
    abs(prev['sn'] - new['sn']) > 0.0001:
    prev = new
    new = getShift(easting - prev['se'], northing - prev['sn'])
  E = easting - new['se']
  N = northing - new['sn']

  a = 6378137.000
  b =  6356752.3141
  e2 = (a**2-b**2)/a**2
  F0 = 0.9996012717
  lat0 = 49 * (math.pi/180)
  lon0 = -2 * (math.pi/180)
  E0 = 400000
  N0 = -100000

  lat_dash = ((N-N0)/(a*F0)) + lat0
  n = (a-b)/(a+b)
  M = b*F0*((1+n+5.0/4.0*n**2+5.0/4.0*n**3)*(lat_dash-lat0)-(3*n+3*n**2+21.0 \
    /8.0*n**3)*math.sin(lat_dash-lat0)*math.cos(lat_dash+lat0)+(15.0/8.0*n** \
    2+15.0/8.0*n**3)*math.sin(2*(lat_dash-lat0))*math.cos(2*(lat_dash+lat0)) \
    -35.0/24.0*n**3*math.sin(3*(lat_dash-lat0))*math.cos(3*(lat_dash+lat0)))
  while abs(N - N0 - M) >= 0.01:
    lat_dash = ((N-N0-M)/(a*F0))+lat_dash
    M = b*F0*((1+n+5.0/4.0*n**2+5.0/4.0*n**3)*(lat_dash-lat0)-(3*n+3*n**2+ \
      21.0/8.0*n**3)*math.sin(lat_dash-lat0)*math.cos(lat_dash+lat0)+(15.0/ \
      8.0*n**2+15.0/8.0*n**3)*math.sin(2*(lat_dash-lat0))*math.cos(2*( \
      lat_dash+lat0))-35.0/24.0*n**3*math.sin(3*(lat_dash-lat0))*math.cos(3* \
      (lat_dash+lat0)))
  
  v = a*F0*(1-(e2*(math.sin(lat_dash)**2)))**-0.5
  p = a*F0*(1-e2)*(1-(e2*(math.sin(lat_dash)**2)))**-1.5
  n2 = v/p-1
  VII = math.tan(lat_dash)/(2*p*v)
  VIII = (math.tan(lat_dash)/(24*p*v**3))*(5+3*(math.tan(lat_dash)**2)+n2-9* \
    (math.tan(lat_dash)**2)*n2)
  IX = (math.tan(lat_dash)/(720*p*v**5))*(61+90*(math.tan(lat_dash)**2)+45*( \
    math.tan(lat_dash)**4))
  X = sec(lat_dash)/v
  XI = (sec(lat_dash)/(6*v**3))*((v/p)+2*(math.tan(lat_dash)**2))
  XII = sec(lat_dash)/(120*v**5)*(5+28*(math.tan(lat_dash)**2)+24*(math.tan( \
    lat_dash)**4))
  XIIA = sec(lat_dash)/(5040*v**7)*(61+662*(math.tan(lat_dash)**2)+1320*( \
    math.tan(lat_dash)**4)+720*(math.tan(lat_dash)**6))

  lat = lat_dash - VII*(E-E0)**2 + VIII*(E-E0)**4 - IX*(E-E0)**6
  lng = lon0 + X*(E-E0) - XI*(E-E0)**3 + XII*(E-E0)**5 - XIIA*(E-E0)**7

  lat = float('%.8f'%(lat * 180 / math.pi))
  lng = float('%.8f'%(lng * 180 / math.pi))

  return {'EN':[easting, northing], 'LL':[lat, lng]}

# Test of conversion from ETRS89 to OSGB36 and back again:

OSGB36 = ETRS89ToOSGB36(51.292798, -0.793407) 
print(OSGB36)

ETRS89 = OSGB36ToETRS89(484228.591, 155542.948)
print(ETRS89)

