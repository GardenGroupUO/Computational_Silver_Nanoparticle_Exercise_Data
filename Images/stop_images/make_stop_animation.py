# import Image, ImageDraw

# im = Image.open("Two_Dalmatians.jpg")

# draw = ImageDraw.Draw(im)

# # Locate the "moon" in the upper-left region of the image
# xy=[x/4 for x in im.size+im.size]

# # Bounding-box is 40x40, so radius of interior circle is 20
# xy=[xy[0]-20, xy[1]-20, xy[2]+20, xy[3]+20]

# # Fill a chord that starts at 45 degrees and ends at 225 degrees.
# draw.chord(xy, 45, 45+180, outline="white", fill="white")

# del draw

# # save to a different file
# with open("Two_Dalmatians_Plus_Moon.png", "wb") as fp:
#     im.save(fp, "PNG")



from PIL import Image, ImageDraw
size = 500
size = (size, size)
centre = (size[0]/2,size[1]/2)

def add_square(draw,centre,length,color):
    xx,yy = centre
    radius = length/2.0
    position = (xx-radius, yy-radius, xx+radius, yy+radius)
    draw.rectangle(position, fill=color)

def add_circle(draw,centre,radius,color):
    xx,yy = centre
    position = (xx-radius, yy-radius, xx+radius, yy+radius)
    draw.ellipse(position, fill=color)

def add_semicircle(draw,centre,radius,angle,width,color):
    xx,yy = centre
    position = (xx-radius, yy-radius, xx+radius, yy+radius)
    draw.arc(position, start=angle, end=angle+180, fill=color, width=width)

color_circle = (66,66,66)

images = []
for angle in range(0,360,10):
    img = Image.new('RGBA', size, (255, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    add_circle(draw,centre=centre,radius=int(size[0]*42/130),color=color_circle)
    add_square(draw,centre=centre,length=int(size[0]*24/130),color=(224,224,224))
    width=int(size[0]*13/130)
    add_semicircle(draw,centre=centre,radius=int(size[0]*60/130),angle=angle,width=width,color=color_circle)
    images.append(img)

images[0].save('stopsvg.gif',save_all=True, append_images=images[1:], optimize=False,disposal=2 , duration=35, loop=0, transparency=0)