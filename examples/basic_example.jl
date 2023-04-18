using OrthoView

# display an interactive OrthoView window with random data 
ortho(rand(30,40,50))

using TestImages

obj = Float32.(testimage("simple_3d_ball"));
ortho(obj)

# not working:
# obj = testimage("simple_3d_ball");
# ortho(obj)

obj = testimage("");
