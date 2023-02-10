# a,b,c,d are 4 points in 2D plane

area_doubled = a * b + b * c + c * d + d * a #(This is for rectangle only)
border_points = math.gcd(a, b) + math.gcd(b, c) + math.gcd(c, d) + math.gcd(d, a) #Int
inner_points = 1 + (area_doubled - border_points) / 2
