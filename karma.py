from manim import *
from math import *
from manim.opengl import *

"""
width=int(1080)
height=int(1920)
config.frame_size = [width, height]
"""

def get_2D_lattice(
        point_count_h, point_count_v, 
        shape = Cube,
        side_length = 0.2, 
        i_ah = 1, i_av = 1,
        real_color = BLUE,
        recip_color = RED,

        automatic = True
):
    lattice_1D = VGroup()
    ah_lattice = ValueTracker(i_ah) #Lattice constant horizontal
    av_lattice = ValueTracker(i_av) #Lattice constant vertical

    lattice_center = [
        ValueTracker(0),
        ValueTracker(0),
        ValueTracker(0),
    ]
    for i in range(point_count_h):
        start = -point_count_h/2 * ah_lattice.get_value()
        point = shape(
            side_length = side_length
        ).set_color(real_color).move_to(RIGHT*(start + ah_lattice.get_value()*i))
        lattice_1D.add(point)
    lattice_1D = lattice_1D.move_to([0,0,0])

    lattice_2D = VGroup()

    for i in range(point_count_v):
        start = -point_count_v/2 * av_lattice.get_value()
        lattice_2D.add(lattice_1D.copy().move_to(DOWN*(start + av_lattice.get_value()*i)))
    lattice_2D.move_to([0,0,0
    ])

    if automatic == False:
        return lattice_2D

    for i in range(point_count_v):
        for j in range(point_count_h):
            lattice_2D[i][j].add_updater(
                lambda mob, i=i, j=j: mob.move_to([
                    lattice_center[0].get_value() + ah_lattice.get_value()*(j - int(point_count_h/2)),
                    lattice_center[1].get_value() + av_lattice.get_value()*(i - int(point_count_v/2)),
                    lattice_center[2].get_value()
                ])
            )

    #RECIPROCAL LATTICE
    recip_1D = VGroup()
    recip_center = [
        ValueTracker(0),
        ValueTracker(0),
        ValueTracker(0),
    ]
    for i in range(point_count_h):
        start = -point_count_h/2 * 1/ah_lattice.get_value()
        point = shape(
            side_length = side_length
        ).set_color(recip_color).move_to(RIGHT*(start + 1/ah_lattice.get_value()*i))
        recip_1D.add(point)
    recip_1D = recip_1D.move_to([0,0,0])

    recip_2D = VGroup()

    for i in range(point_count_v):
        start = -point_count_v/2 * 1/av_lattice.get_value()
        recip_2D.add(recip_1D.copy().move_to(DOWN*(start + 1/av_lattice.get_value()*i)))
    recip_2D.move_to([0,0,0])

    for i in range(point_count_v):
        for j in range(point_count_h):
            recip_2D[i][j].add_updater(
                lambda mob, i=i, j=j: mob.move_to([
                    recip_center[0].get_value() + 1/ah_lattice.get_value()*(j - int(point_count_h/2)),
                    recip_center[1].get_value() + 1/av_lattice.get_value()*(i - int(point_count_v/2)),
                    recip_center[2].get_value()
                ])
            )

    return {
        'real': [lattice_2D, lattice_center],
        'reciprocal': [recip_2D, recip_center],
        'constants': [ah_lattice, av_lattice]
    }

def get_optimal_2D_lattice(point_count_h, point_count_v, r = 0.1, horiz = 1, vert = 0.5):
    lattice_1D = VGroup()

    lattice_center = [
        ValueTracker(0),
        ValueTracker(0),
        ValueTracker(0),
    ]
    for i in range(point_count_h):
        start = -point_count_h/2 * horiz
        point = Dot(
            radius = r
        ).set_color(BLUE).move_to(RIGHT*(start + horiz*i))
        lattice_1D.add(point)
    lattice_1D = lattice_1D.move_to([0,0,0])

    lattice_2D = VGroup()

    shifted = False
    for i in range(point_count_v):
        start = -point_count_v/2 * vert

        if shifted == False:
            lattice_2D.add(lattice_1D.copy().move_to(DOWN*(start + vert*i)))
            shifted = True
        else:
            lattice_2D.add(lattice_1D.copy().move_to(DOWN*(start + vert*i) + RIGHT*horiz/2).set_color(RED)[:len(lattice_1D)-1])
            shifted = False

    lattice_2D.move_to([0,0,0])

    return lattice_2D

def get_eighth_circle(radius = 1, color = BLUE, resolution = 4, side_opacity = 1):
    r = radius

    # Define the parametric surface for 1/8 of the sphere
    surface = Surface(
        lambda u, v: np.array([
            r * np.sin(u) * np.cos(v),
            r * np.sin(u) * np.sin(v),
            r * np.cos(u)
        ]),
        u_range=[0, PI / 2],
        v_range=[0, PI / 2],
        resolution = resolution
    ).set_color(color)

    sides = VGroup(Sector(
                inner_radius = 0,
                outer_radius = radius,     
                start_angle=0,            
                angle=PI / 2,             
                color=color
            ).set_opacity(side_opacity))
    
    sides.add(sides[0].copy().rotate(angle=PI/2, axis=DOWN, about_point = [0,0,0]))
    sides.add(sides[0].rotate(angle=PI/2, axis=RIGHT, about_point = [0,0,0]))

    return VGroup(sides, surface)

def get_eighth_circle_GL(radius = 1, color = BLUE, side_opacity = 1):
    r = radius

    # Define the parametric surface for 1/8 of the sphere
    surface = OpenGLSurface(
        lambda u, v: np.array([
            r * np.sin(u) * np.cos(v),
            r * np.sin(u) * np.sin(v),
            r * np.cos(u)
        ]),
        u_range=[0, PI / 2],
        v_range=[0, PI / 2],
    ).set_color(color)

    sides = VGroup(Sector(
                inner_radius = 0,
                outer_radius = radius,     
                start_angle=0,            
                angle=PI / 2,             
                color=color
            ).set_opacity(side_opacity))
    
    sides.add(sides[0].copy().rotate(angle=PI/2, axis=DOWN, about_point = [0,0,0]))
    sides.add(sides[0].rotate(angle=PI/2, axis=RIGHT, about_point = [0,0,0]))

    return Group(sides, surface)

class CircleSections(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi = 30*DEGREES
        )

        r = 2

        # Define the parametric surface for 1/8 of the sphere
        surface = OpenGLSurface(
            lambda u, v: np.array([
                r * np.sin(u) * np.cos(v),
                r * np.sin(u) * np.sin(v),
                r * np.cos(u)
            ]),
            u_range=[0, PI / 2],
            v_range=[0, PI / 2],
        ).set_color(BLUE)

        sides = VGroup(Sector(
                    inner_radius = 0,
                    outer_radius = 2,     
                    start_angle=0,            
                    angle=PI / 2,             
                    color=BLUE
                ))
        
        sides.add(sides[0].copy().rotate(angle=PI/2, axis=DOWN, about_point = [0,0,0]))
        sides.add(sides[0].rotate(angle=PI/2, axis=RIGHT, about_point = [0,0,0]))

        # Add the surface to the scene
        self.add(surface)
        self.add(sides)

        self.play(Rotate(Group(surface, sides), angle = PI/4, axis = [0.2,0.8,0]), run_time = 4)

        self.play(self.camera.animate.set_euler_angles(phi = -120*DEGREES))


class KarmaApproachGL(ThreeDScene):
    def construct(self):

        self.set_camera_orientation(
            phi = -20*DEGREES,
            zoom = 1.2
        )
        BCC = VGroup(Cube(
            side_length = 2.5
        ).set_fill(opacity = 0))

        
        corners = [
            [1.25, 1.25, 0],
            [1.25, -1.25, 0],
            [-1.25, -1.25, 0],
            [-1.25, 1.25, 0]
        ]
        for i in range(len(corners)):
            corner = np.array(corners[i])
            corners.append(list(corner - np.array([0,0,2.5])))
        

        for corner in corners:
            BCC.add(
                Sphere(radius=0.1, resolution=15).move_to(corner)
            ).set_color(BLUE)
        BCC.add(Sphere(radius=0.1, resolution=15).set_color(RED).move_to([0,0,-1.25]))

        self.play(Create(BCC))

        self.wait()

        self.play(
            BCC[0].move_to(BCC[-1].get_center()).animate.set_stroke(WHITE, width=2, opacity = 0.5)
        )
        self.wait()
        self.play(BCC.animate.move_to([0,0,0]))
        self.wait()


        self.play(
            *[
                sphere.animate.scale(10) for sphere in BCC[1:]
            ]
        )

        self.wait()

        self.play(Rotate(BCC, angle = 2*PI, axis=UP), run_time = 10, rate_func=linear)
        self.wait()

        unit_cell = VGroup()
        for corner in corners:
            unit_cell.add(
                get_eighth_circle(color=GOLD, resolution = 15, side_opacity = 0.5).shift(corner)
            ).set_color(GOLD)

        for i in range(int(len(unit_cell)/2)):
            unit_cell[i].rotate(
                angle = PI/2 * (2 - i),
                axis = OUT,
                about_point = corners[i]
            ).rotate(
                angle = copysign(PI/2, unit_cell[i].get_center()[1]),
                axis = RIGHT,
                about_point = corners[i]
            )
            unit_cell[i+4].rotate(
                angle = PI/2 * (2 - i),
                axis = OUT,
                about_point = corners[i+4]
            )
        unit_cell.add(Sphere(radius=1, resolution=15).set_color(GOLD).move_to([0,0,-1.25]))

        unit_cell.move_to([0,0,0])
        self.play(FadeIn(unit_cell), BCC.animate.set_opacity(0.2))
        self.wait()
        self.play(Rotate(VGroup(BCC, unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(FadeOut(BCC))
        self.wait()
        self.play(Rotate(VGroup(unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(*[
            mob.animate.shift(
                -np.array(corners[i]) + IN*1.25
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
        ])
        self.wait()
        self.play(Rotate(VGroup(unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(*[
            mob.animate.shift(
                np.array(corners[i]) + OUT*1.25
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
        ])
        self.wait()

class KarmaApproach(ThreeDScene):
    def construct(self):

        self.set_camera_orientation(
            phi = -20*DEGREES,
            zoom = 1.2
        )
        BCC = VGroup(Cube(
            side_length = 2.5
        ).set_fill(opacity = 0))

        
        corners = [
            [1.25, 1.25, 0],
            [1.25, -1.25, 0],
            [-1.25, -1.25, 0],
            [-1.25, 1.25, 0]
        ]
        for i in range(len(corners)):
            corner = np.array(corners[i])
            corners.append(list(corner - np.array([0,0,2.5])))
        

        for corner in corners:
            BCC.add(
                Sphere(radius=0.1, resolution=15).move_to(corner)
            ).set_color(BLUE)
        BCC.add(Sphere(radius=0.1, resolution=15).set_color(RED).move_to([0,0,-1.25]))

        self.play(Create(BCC))

        self.wait()

        self.play(
            BCC[0].move_to(BCC[-1].get_center()).animate.set_stroke(WHITE, width=2, opacity = 0.5)
        )
        self.wait()
        self.play(BCC.animate.move_to([0,0,0]))
        self.wait()


        self.play(
            *[
                sphere.animate.scale(10) for sphere in BCC[1:]
            ]
        )

        self.wait()

        self.play(Rotate(BCC, angle = 2*PI, axis=UP), run_time = 10, rate_func=linear)
        self.wait()

        unit_cell = VGroup()
        for corner in corners:
            unit_cell.add(
                get_eighth_circle(color=GOLD, resolution = 15, side_opacity = 0.5).shift(corner)
            ).set_color(GOLD)

        for i in range(int(len(unit_cell)/2)):
            unit_cell[i].rotate(
                angle = PI/2 * (2 - i),
                axis = OUT,
                about_point = corners[i]
            ).rotate(
                angle = copysign(PI/2, unit_cell[i].get_center()[1]),
                axis = RIGHT,
                about_point = corners[i]
            )
            unit_cell[i+4].rotate(
                angle = PI/2 * (2 - i),
                axis = OUT,
                about_point = corners[i+4]
            )
        unit_cell.add(Sphere(radius=1, resolution=15).set_color(GOLD).move_to([0,0,-1.25]))

        unit_cell.move_to([0,0,0])
        self.play(FadeIn(unit_cell), BCC.animate.set_opacity(0.2))
        self.wait()
        self.play(Rotate(VGroup(BCC, unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(FadeOut(BCC))
        self.wait()
        self.play(Rotate(VGroup(unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(*[
            mob.animate.shift(
                -np.array(corners[i]) + IN*1.25
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
        ])
        self.wait()
        self.play(Rotate(VGroup(unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(*[
            mob.animate.shift(
                np.array(corners[i]) + OUT*1.25
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
        ])
        self.wait()

class TwoDLattice(Scene):
    def construct(self):

        sides = Square(side_length = 2.5).set_stroke(opacity=0.5)
        side_length = MathTex("a").next_to(sides, LEFT)
        side_length = VGroup(BackgroundRectangle(side_length), side_length)

        corners = [
            [1.25, 1.25, 0],
            [1.25, -1.25, 0],
            [-1.25, 1.25, 0],
            [-1.25, -1.25, 0]
        ]

        lattice = VGroup()
        for corner in corners:
            lattice.add(
                Dot().move_to(corner).set_color(BLUE)
            )
        lattice.add(Dot().set_color(RED))

        self.play(Create(lattice))
        self.wait()

        self.play(Create(VGroup(sides, side_length)))
        self.wait()

        whole_area = VGroup(MathTex("a^2"))
        whole_area = VGroup(BackgroundRectangle(whole_area), whole_area)
        self.play(Write(whole_area))
        self.wait()
        self.play(FadeOut(whole_area))
        self.wait()

        connecting_lines = VGroup()
        for corner in corners:
            connecting_lines.add(
                DashedLine(
                    start = lattice[-1].get_center(),
                    end = corner
                ).set_color(PURPLE)
            )

        self.play(Create(connecting_lines))
        self.wait()

        bisecting_lines = VGroup()
        for line in connecting_lines:
            midpoint = (line.get_end() + line.get_start())/2
            direction = line.get_end() + line.get_start()

            perpendicular_direction = np.array([-direction[1], direction[0], 0])
            perpendicular_line_length = 1.5
            p1_perpendicular = midpoint + perpendicular_direction * perpendicular_line_length / 2
            p2_perpendicular = midpoint - perpendicular_direction * perpendicular_line_length / 2

            bisecting_line = Line(
                start = p1_perpendicular,
                end = p2_perpendicular,
            ).set_opacity(0.4)

            bisecting_lines.add(bisecting_line)

        self.play(Create(bisecting_lines))
        self.wait()

        fill_area = Square(
            side_length = 1.25*np.sqrt(2)
        ).set_stroke(opacity = 0).set_fill(opacity = 0.2, color = LIGHTER_GRAY).rotate(45*DEGREES)

        self.play(Create(fill_area))
        self.wait()

        domain_sl = MathTex("\\frac{a}{\\sqrt{2}}").next_to(bisecting_lines[2].get_center(), UL, buff = -0.15).scale(0.6)
        domain_sl = VGroup(BackgroundRectangle(domain_sl).set_opacity(0.4), domain_sl)
        self.play(Write(domain_sl))
        self.wait()

        twoD_area = MathTex("\\frac{a^2}{2}")
        twoD_area = VGroup(BackgroundRectangle(twoD_area, buff = 0.1).set_opacity(0.3), twoD_area)
        self.play(Write(twoD_area))
        self.wait()

        #Thinking outside the box
        self.play(FadeOut(domain_details := VGroup(bisecting_lines, connecting_lines, side_length, domain_sl, twoD_area, fill_area)))
        self.wait()
        full_lattice = get_optimal_2D_lattice(10, 15, horiz = 2.5, vert = 1.25)

        self.play(DrawBorderThenFill(full_lattice))
        self.wait()

        center_corner_rects = VGroup(
            SurroundingRectangle(lattice[-1])
        )
        self.play(FadeIn(center_corner_rects))
        self.wait()
        new_sides = VGroup(sides.copy().shift(DOWN*1.25 + RIGHT*1.25))
        self.play(
            ReplacementTransform(sides.copy(), new_sides)
        )
        self.wait()

        center_corner_rects.add(SurroundingRectangle(lattice[1]))
        self.play(FadeIn(center_corner_rects[1]))
        self.wait()

        for point in lattice:
            center_corner_rects.add(SurroundingRectangle(point))
        for corner in corners:
            new_sides.add(sides.copy().move_to(corner))

        self.play(FadeIn(VGroup(center_corner_rects[2:], new_sides[1:])))
        new_sides.add(sides)

        self.wait()

        other_connections = VGroup()
        for i in range(-10, 10, 1):
            for j in range(-10, 10, 1):
                other_connections.add(
                    sides.copy().move_to(
                        [i*1.25, j*1.25, 0]
                    ).set_stroke(opacity=0.2)
                )
        self.play(
            FadeIn(other_connections),
            FadeOut(new_sides)
        )
        self.wait()

        #EXPANDING THE DOTS (TRANSITION INTO NEXT SCENE)

        self.play(
            FadeOut(other_connections, center_corner_rects))
        self.wait()

        self.play(FadeIn(domain_details))
        self.wait()
        self.play(domain_details.animate.shift(UP*1.25 + RIGHT*1.25))
        self.wait()
        self.play(domain_details.animate.shift(DOWN*2.5 + LEFT*2.5))
        self.wait()
        self.play(domain_details.animate.shift(UP*1.25 + LEFT*1.25))
        self.wait()
        self.play(domain_details.animate.shift(RIGHT*2.5))
        self.wait()
        self.play(FadeOut(domain_details))
        self.wait()
        

        self.play(FadeOut(full_lattice))
        self.wait()
        self.play(FadeIn(sides))
        self.wait()
        self.play(
            *[
                mob.animate.scale(10) for mob in lattice
            ]
        )
        self.wait()

        unit_cell = VGroup()
        for i, corner in enumerate(corners):
            unit_cell.add(
                Sector(
                    inner_radius = 0,
                    outer_radius = np.array(lattice[0].get_top() - lattice[0].get_center())[1],                 
                    start_angle=0,            
                    angle=PI / 2,             
                    color=GOLD            
                ).shift(corner)
            )
        
        unit_cell[0].rotate(PI, about_point = corners[0])
        unit_cell[1].rotate(PI/2, about_point = corners[1])
        unit_cell[2].rotate(-PI/2, about_point = corners[2])
        unit_cell.add(Dot().scale(10).set_color(GOLD))
        temp_cell = unit_cell.copy()

        self.play(FadeIn(unit_cell))
        self.wait()

        self.play(
            *[
                Rotate(mob, angle = PI, about_point = corners[i]) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
            ]
        )
        self.play(
            *[
                mob.animate.shift(-np.array(corners[i])) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
            ]
        )

        self.wait()
        self.play(
            *[
                Rotate(mob, angle = PI, about_point = np.array(corners[i])/2) for i, mob in enumerate(unit_cell[:len(unit_cell)-1])
            ]
        )
        self.wait()
        self.play(FadeIn(domain_details))
        self.wait()
        self.play(domain_details.animate.shift(DOWN*1.25 + LEFT*1.25))
        for i in range(len(corners)):
            self.play(
                domain_details.animate.shift(np.array(corners[i]) - np.array(corners[i-1]))
            )
        self.wait()

        self.play(domain_details.animate.shift(corners[0]))
        self.wait()
        self.play(VGroup(lattice, unit_cell, domain_details, sides).animate.shift(DOWN*1))
        self.wait()
        crux = VGroup(
            Tex("The ratio of areas associated with one point"),
            Tex("to another lattice point must be exactly 1:1")
        ).arrange(DOWN).shift(UP*2)
        self.play(Write(crux), run_time = 5)
        self.wait()
        self.play(VGroup(lattice, unit_cell, domain_details, sides).animate.shift(RIGHT*2.5))
        self.wait()

        temp = Dot(color=GOLD).scale(10)
        final_explanation = VGroup(
            Tex("2 total points"),
            VGroup(temp,
                   VGroup(*[mob.rotate(angle = PI, about_point = np.array(corners[i])/2) for i, mob in enumerate(temp_cell[:len(unit_cell)-1])])).arrange(RIGHT).move_to(temp.get_center()).scale(0.707),
            Tex("$A_{tot} = A_{close} + A_{far} = a^2$"),
            MathTex("A_{close} = A_{far} = \\frac{a^2}{2}")
        ).arrange(DOWN).move_to([-3, -1, 0])
        for i, item in enumerate(final_explanation):
            self.play(Write(item))
            if i == 1:
                self.play(
                    *[
                        mob.animate.shift(np.array(corners[i])/4) for i, mob in enumerate(item[1])
                    ]
                )
                self.wait()
                self.play(
                    *[
                        mob.animate.shift(-np.array(corners[i])/4) for i, mob in enumerate(item[1])
                    ]
                )
            self.wait()
        self.wait()


        self.play(
            FadeOut(VGroup(final_explanation, unit_cell, domain_details), crux),
        )
        self.play(VGroup(lattice, sides).animate.move_to(ORIGIN))
        self.play(*[
            mob.animate.scale(0.1) for mob in lattice
        ])
        self.wait()

        #Hexagonal Lattice
        hex_r = 1
        hex_corners = [[hex_r * np.cos(theta), hex_r * np.sin(theta), 0] for theta in np.linspace(0, 2 * np.pi, 7)]  
        true_r = hex_r * np.sqrt(3)/2

        hexagon = Polygon(*hex_corners)

        hex_connections = VGroup()

        dx = np.sqrt(3) * true_r
        dy = 2*true_r  

        for i in range(-5, 5, 1):
            for j in range(-5, 5, 1):
                x = j*dx
                y = i*dy

                if j % 2 == 1:
                    y += true_r
                hex_connections.add(
                    hexagon.copy().move_to(
                        [x, y, 0]
                    )
                )
        hex_connections.set_stroke(opacity=0.3)

        hex_lattice = VGroup()
        for hex in hex_connections:
            for angle in range(6):
                hex_lattice.add(
                    Dot().move_to(
                        hex.get_center() + RIGHT*hex_r
                    ).rotate(angle = angle*PI/3, about_point = hex.get_center()).set_color(BLUE)
                )
        
        center_point = hex_lattice[50]
        for point in hex_lattice:
            if np.linalg.norm(center_point.get_center()) > np.linalg.norm(point.get_center()):
                center_point = point

        VGroup(hex_lattice, hex_connections).shift(-center_point.get_center())

        self.play(ReplacementTransform(lattice, hex_lattice), FadeOut(sides))
        self.wait()
        self.play(Create(hex_connections))
        self.wait()

        self.play(VGroup(hex_connections, hex_lattice).animate.scale(3, about_point = center_point.get_center()))
        self.wait()

        self.play(
            mob.animate.scale(5) for mob in hex_lattice
        )
        self.wait()
        self.play(
            mob.animate.scale(1/5) for mob in hex_lattice
        )
        self.wait()

        hex_bisectors = VGroup()
        for i in range(3):
            connection = DashedLine(
                    start = center_point.get_center(),
                    end = hex_r*3 * np.array([np.cos(2*PI/3*i), np.sin(2*PI/3*i), 0])
            ).set_opacity(0.5).set_color(PURPLE)

            hex_bisectors.add(
                connection.copy().rotate(PI/2).scale(2)
            )
        self.play(Create(hex_bisectors))
        self.wait()

        hex_cell = hexagon.copy().scale(3).set_stroke(opacity = 0.4).set_color(GOLD).move_to(center_point.get_center())
        self.play(FadeIn(hex_cell))
        self.wait()

class Testing(ThreeDScene):
    def construct(self):
        l = np.sqrt((0.3*2.92)**2 + (2.5*np.sqrt(3)/4)**2)
        a = np.sqrt((0.25*0.3*2.92*np.sqrt(3))**2 + l**2)
        b = 2.5*np.sqrt(2)-2*a
        c = b*np.sqrt(2)
        h = 2.5 - 2*c

        wedge_vertices = [
            [-1.25,-1.25,0],
            [-0.625, 0, -h],
            [-1.25,1.25,0],
            [-1.25+h, 0, 0]
        ]
        wedge_faces = [
            [0, 1, 2],
            [0, 1, 3],
            [1, 2, 3],
            [0, 2, 3]
        ]
        wedge = Polyhedron(
            wedge_vertices,
            wedge_faces
        )

        wedges = VGroup(
            *[
                wedge.copy().rotate(about_point=[0,0,0], angle = PI/2*i) for i in range(4)
            ]
        )

        flipped = wedges.copy().rotate(
            angle = PI,
            about_point = [0,0, -1.25],
            axis = UP
        )
        self.play(Create(VGroup(wedges, flipped)))

        vert_wedge_vertices = [
            [-1.25,-1.25,0],
            [-1.25+h,-1.25,-1.25],
            [-1.25,-1.25,-2.5],
            [-1.25,-1.25+h,-1.25]
        ]
        vert_wedge = Polyhedron(
            vert_wedge_vertices,
            wedge_faces
        )
        vert_wedges = VGroup(
            *[
                vert_wedge.copy().rotate(about_point=[0,0,-1.25], angle = PI/2*i, axis = OUT) for i in range(4)
            ]
        )

        wedges = VGroup(
            *[mob for mob in wedges],
            *[mob for mob in flipped],
            *[mob for mob in vert_wedges]
        )

        self.wait()
        self.play(Rotate(wedges, axis = [0.4, 0.8, 0], angle = 2*PI), run_time = 10)
        self.wait()

class BCCLattice(ThreeDScene):
    def construct(self):

        self.set_camera_orientation(
            phi = -20*DEGREES,
            zoom = 2
        )
        BCC = VGroup(Cube(
            side_length = 2.5
        ).set_fill(opacity = 0))

        
        corners = [
            [1.25, 1.25, 0],
            [1.25, -1.25, 0],
            [-1.25, 1.25, 0],
            [-1.25, -1.25, 0]
        ]
        for i in range(len(corners)):
            corner = np.array(corners[i])
            corners.append(list(corner - np.array([0,0,2.5])))
        

        for corner in corners:
            BCC.add(
                Sphere(radius=0.1, resolution=6).move_to(corner)
            )
        BCC.add(Sphere(radius=0.1, resolution=6).set_color(RED).move_to([0,0,-1.25]))

        self.play(Create(BCC))

        self.wait()

        self.play(
            BCC[0].move_to(BCC[-1].get_center()).animate.set_stroke(WHITE, width=2, opacity = 0.5)
        )
        self.wait()

        connecting_lines = VGroup()
        for corner in BCC[1:len(BCC)-1]:
            connecting_lines.add(
                Line(
                    start = corner.get_center(),
                    end = BCC[-1].get_center()
                ).set_opacity(0.5).set_color(PURPLE)
            )

        self.play(Create(connecting_lines))
        self.wait()

        #Creating Planes
        bisecting_planes = VGroup()
        bisecting_hexagons = VGroup()
        bisecting_dict = {}
        for i, line in enumerate(connecting_lines):
            midpoint = (line.get_end() - line.get_start())/2
            normal_vector = line.get_end() - line.get_start()
            normal_vector = normal_vector/np.linalg.norm(normal_vector)
            current_normal = OUT
            rotation_angle = angle_between_vectors(current_normal, normal_vector)
            rotation_axis = np.cross(current_normal, normal_vector)
            
            radius = 0.3 * 2.92
            hex_corners = [
                [radius * np.cos(theta), radius * np.sin(theta), 0]
                for theta in np.linspace(0, 2 * np.pi, 7)  # 7 points to close the loop
            ]

            start_corner = hex_corners[0]
            start_to_corner =  -np.array(start_corner)
            center_to_corner = -np.array([line.get_start()[0], line.get_start()[1], 0])
            initial_angle = angle_between_vectors(start_to_corner, center_to_corner)

            hexagon = Polygon(*hex_corners).rotate(initial_angle)   
            hexagon.rotate(rotation_angle, axis=rotation_axis, about_point=ORIGIN).move_to(midpoint).set_opacity(0.3).set_color(GREEN)

            square = Square(side_length=0.4).set_color(GREEN).rotate(rotation_angle, axis=rotation_axis, about_point=ORIGIN).move_to(midpoint).set_opacity(0.5).set_color(GREEN)

            if hexagon.get_center()[1] < 0:
                hexagon.rotate(angle = PI/6, axis = normal_vector)

            bisecting_hexagons.add(hexagon)
            bisecting_planes.add(square)
            key = f'{i}'
            bisecting_dict[key] = hexagon

        self.play(
            Create(bisecting_planes.move_to(BCC[-1].get_center())),
            run_time = 6
        )

        self.set_camera_orientation(
            frame_center = [0,0,0]
        )
        self.play(
            VGroup(BCC, connecting_lines, bisecting_planes).animate.move_to([0,0,0])
        )
        self.wait()

        self.play(
            Rotate(VGroup(BCC, connecting_lines, bisecting_planes), axis = UP, angle = 3*PI/3),
            rate_func = linear,
            run_time = 6
        )

        self.wait()

        for i in range(len(bisecting_planes)):
            closest_hex = None
            closest_distance = 10
            for hexagon in bisecting_hexagons:
                if np.linalg.norm(hexagon.get_center() - bisecting_planes[i].get_center()) < closest_distance:
                    closest_hex = hexagon
                    closest_distance = np.linalg.norm(hexagon.get_center() - bisecting_planes[i].get_center())
            self.play(ReplacementTransform(bisecting_planes[i], closest_hex))

        self.wait()

        self.play(
            Rotate(VGroup(BCC, connecting_lines, bisecting_hexagons), axis = UP, angle = 3*PI/3),
            rate_func = linear,
            run_time = 20
        )

        self.wait()

        #ADDING NEGATIVE SPACE
        hex_vertices = []
        for i in range(6):
            hex_vertices.append(
                0.3*2.92*np.array(
                    [np.cos(PI/3*i), np.sin(PI/3*i), 0]
                )
            )   
        hex_vertices.append(np.array([0,0,2.5*np.sqrt(3)/4]))
        hex_faces = []
        for i in range(len(hex_vertices)-1):
            if i == 0:
                hex_faces.append(
                    [len(hex_vertices)-2, i, len(hex_vertices)-1]
                )
            else:   
                hex_faces.append(
                    [i-1, i, len(hex_vertices)-1]
                )

        hex_pyramid = Polyhedron(
            hex_vertices,
            hex_faces
        )
        for i in range(len(hex_pyramid.faces)):
            if i%2 == 0:
                hex_pyramid.faces[i].set_color(PURPLE)

        hex_pyramids = VGroup()
        for i, line in enumerate(connecting_lines):
            midpoint = (line.get_end() - line.get_start())/2
            normal_vector = line.get_end() - line.get_start()
            normal_vector = normal_vector/np.linalg.norm(normal_vector)
            current_normal = OUT
            rotation_angle = angle_between_vectors(current_normal, normal_vector)
            rotation_axis = np.cross(current_normal, normal_vector)

            start_corner = [0.3*2.92, 0, 0]
            start_to_corner =  -np.array(start_corner)
            center_to_corner = -np.array([line.get_start()[0], line.get_start()[1], 0])
            initial_angle = angle_between_vectors(start_to_corner, center_to_corner)

            temp_pyramid = hex_pyramid.copy().rotate(initial_angle, axis=OUT, about_point=[0,0,0])
            temp_pyramid.rotate(rotation_angle, axis=rotation_axis, about_point=ORIGIN).move_to(midpoint)

            if temp_pyramid.get_center()[1] < 0:
                temp_pyramid.rotate(angle = PI/6, axis = normal_vector)

            hex_pyramids.add(temp_pyramid)

        self.play(Create(hex_pyramids))
        self.wait()
        self.play(FadeOut(bisecting_hexagons))
        self.wait()

        self.play(
            Rotate(VGroup(hex_pyramids, BCC, connecting_lines), axis = UP, angle = 6*PI/3),
            rate_func = linear,
            run_time = 12
        )
        self.wait()

        l = np.sqrt((0.3*2.92)**2 + (2.5*np.sqrt(3)/4)**2)
        a = np.sqrt((0.25*0.3*2.92*np.sqrt(3))**2 + l**2)
        b = 2.5*np.sqrt(2)-2*a
        c = b*np.sqrt(2)
        h = 2.5 - 2*c

        wedge_vertices = [
            [-1.25,-1.25,0],
            [-0.625, 0, -h],
            [-1.25,1.25,0],
            [-1.25+h, 0, 0]
        ]
        wedge_faces = [
            [0, 1, 2],
            [0, 1, 3],
            [1, 2, 3],
            [0, 2, 3]
        ]
        wedge = Polyhedron(
            wedge_vertices,
            wedge_faces
        )

        wedges = VGroup(
            *[
                wedge.copy().rotate(about_point=[0,0,0], angle = PI/2*i) for i in range(4)
            ]
        )

        flipped = wedges.copy().rotate(
            angle = PI,
            about_point = [0,0, -1.25],
            axis = UP
        )

        vert_wedge_vertices = [
            [-1.25,-1.25,0],
            [-1.25+h,-1.25,-1.25],
            [-1.25,-1.25,-2.5],
            [-1.25,-1.25+h,-1.25]
        ]
        vert_wedge = Polyhedron(
            vert_wedge_vertices,
            wedge_faces
        )
        vert_wedges = VGroup(
            *[
                vert_wedge.copy().rotate(about_point=[0,0,-1.25], angle = PI/2*i, axis = OUT) for i in range(4)
            ]
        )

        wedges = VGroup(
            *[mob for mob in wedges],
            *[mob for mob in flipped],
            *[mob for mob in vert_wedges]
        ).set_color(RED).shift(OUT*1.25)

        self.play(Create(wedges))
        self.wait()
        self.play(
            Rotate(VGroup(hex_pyramids, BCC, connecting_lines, wedges), axis = UP, angle = 6*PI/3),
            rate_func = linear,
            run_time = 12
        )
        self.wait()


class NegativeSpace(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=20*DEGREES, zoom = 1.5)

        hex_vertices = []
        for i in range(6):
            hex_vertices.append(
                0.3*2.92*np.array(
                    [np.cos(PI/3*i), np.sin(PI/3*i), 0]
                )
            )   
        hex_vertices.append(np.array([0,0,2.5*np.sqrt(3)/4]))
        hex_faces = []
        for i in range(len(hex_vertices)-1):
            if i == 0:
                hex_faces.append(
                    [len(hex_vertices)-2, i, len(hex_vertices)-1]
                )
            else:   
                hex_faces.append(
                    [i-1, i, len(hex_vertices)-1]
                )

        hex_pyramid = Polyhedron(
            hex_vertices,
            hex_faces
        )

        l = np.sqrt((0.3*2.92)**2 + (2.5*np.sqrt(3)/4)**2)
        a = np.sqrt((0.25*0.3*2.92*np.sqrt(3))**2 + l**2)
        b = 2.5*np.sqrt(2)-2*a
        c = b*np.sqrt(2)
        h = 2.5 - 2*c

        wedge_vertices = [
            [-1.25,-1.25,0],
            [-0.625, 0, -h],
            [-1.25,1.25,0],
            [-1.25+h, 0, 0]
        ]
        wedge_faces = [
            [0, 1, 2],
            [0, 1, 3],
            [1, 2, 3],
            [0, 2, 3]
        ]
        wedge = Polyhedron(
            wedge_vertices,
            wedge_faces
        )

        for i in range(len(hex_pyramid.faces)):
            if i%2 == 0:
                hex_pyramid.faces[i].set_color(PURPLE)
        for i in range(len(wedge.faces)):
            wedge.faces[i].set_color(RED)

        self.play(Create(VGroup(
            hex_pyramid.shift(LEFT*1.5+IN*0.5),
            wedge.shift(RIGHT*2).rotate(PI/1.2 + 4*PI/5, axis=UP)
        )))
        self.wait()

        self.move_camera(
            phi = 80*DEGREES,
            run_time = 5
        )
        self.wait()
        self.play(
            Rotate(VGroup(wedge, hex_pyramid), about_point=[0,0,0],axis=OUT),
            run_time = 5
        )
        self.wait()

class NegativeSpaceCalculations(Scene):
    def construct(self):
        derivation = VGroup(
            Tex("$h$: Pyramid Height, $L$: Octohedron Edge, $a$: Lattice Constant, $l$: Base Edge").scale(0.6),
            VGroup(
                VGroup(
                    Text("Hexagonal Pyramid", font='Arial'),
                    Tex("Volume is the standard formula"),
                    MathTex("V_p = \\frac{\\sqrt{3}}{2}l^2h,", "\\text{ }", "h = \\frac{a\\sqrt{3}}{4},", "l = \\frac{3}{6\\sqrt{2}}a"),
                    MathTex("V_p = \\frac{\\sqrt{3}}{2}\\frac{9}{72}a^2\\frac{\\sqrt{3}}{4}a", "=", "\\frac{3}{64}a^3"),
                    Tex("8 pyramids, so the total volume is $\\frac{3}{8}a^3$"),
                ).arrange(DOWN).scale(0.6),
                VGroup(
                    Text("Wedge", font='Arial'),
                    Tex("Split the wedge into 2 tetrahedrons"),
                    MathTex("V_{tet} = \\frac{1}{6}\\text{(l)(w)(h)}, \\text{h}=\\frac{a}{2}, \\text{l} = \\sqrt{\\text{l}^2+\\text{h}^2-\\frac{a^2}{4}}").scale(0.9),
                    # MathTex("V_{tet} = \\frac{a}{2}\\left(\\frac{9}{72}a^2+\\frac{3}{16}a^2-\\frac{a^2}{4}\\right)^2"),
                    MathTex("V_{tet} = \\frac{1}{6}\\left(\\frac{9}{144}a^3+\\frac{3}{32}a^3-\\frac{a^3}{8}\\right)} = \\frac{a^3}{192}"),
                    Tex("24 tetrahedrons, we have $\\frac{24}{192} = \\frac{1}{8}a^3$")
                ).arrange(DOWN).scale(0.6)
            ).arrange(RIGHT),
            Tex("Put it all together and tada: $\\frac{a^3}{2}$")
        ).arrange(DOWN, buff = 0.1)

        self.play(Write(derivation))
        self.wait()



class KarmaApproachGL(ThreeDScene):
    def construct(self):

        self.set_camera_orientation(
            phi = -20*DEGREES,
        )

        corners = [
            [1.25, 1.25, 0],
            [1.25, -1.25, 0],
            [-1.25, -1.25, 0],
            [-1.25, 1.25, 0]
        ]
        for i in range(len(corners)):
            corner = np.array(corners[i])
            corners.append(list(corner - np.array([0,0,2.5])))

        unit_cell = Group()
        for corner in corners:
            unit_cell.add(
                get_eighth_circle_GL(color=GOLD, side_opacity = 0.5).shift(corner)
            ).set_color(GOLD)

        for i in range(int(len(unit_cell)/2)):
            unit_cell[i].rotate(
                angle = PI/2 * (2 - i),
                axis = OUT,
                about_point = corners[i]
            ).rotate(
                angle = copysign(PI/2, unit_cell[i].get_center()[1]),
                axis = RIGHT,
                about_point = corners[i]
            )
            unit_cell[i+4].rotate(
                angle = PI/2 * (2 - i),
                axis = OUT,
                about_point = corners[i+4]
            )
        # unit_cell.add(Sphere(radius=1, resolution=15).set_color(GOLD).move_to([0,0,-1.25]))

        unit_cell.move_to([0,0,0])
        unit_cell.add(*[
            mob.copy().shift(
                -np.array(corners[i]) + IN*1.25
            ).set_opacity(0) for i, mob in enumerate(unit_cell[:len(unit_cell)])
        ])
        self.play(Create(unit_cell))
        self.wait()
        # self.play(Rotate(unit_cell, axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        # self.play(Rotate(unit_cell, axis= [0.1,0.9,0], about_point=ORIGIN), run_time = 5)
        self.wait()

        """
        self.play(
            *[Rotate(
                mob,
                axis = np.array([-corners[i][1], corners[i][0], corners[i][2]])/2,
                angle = PI,
                about_point = corner
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-8])]
        )"""

        self.play(*[
            mob.animate.shift(
                -np.array(corners[i]) + IN*1.25
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-8])
        ])
        self.wait()
        # self.play(Rotate(Group(unit_cell), axis= [0.1,0.9,0]), run_time = 5)
        self.wait()
        self.play(*[
            mob.animate.shift(
                np.array(corners[i]) + OUT*1.25
            ) for i, mob in enumerate(unit_cell[:len(unit_cell)-8])
        ])
        self.wait()

        self.interactive_embed()
