from manim import *
import numpy as np

def fade_all(self):
    self.play(
        *[FadeOut(mob) for mob in self.mobjects]
    )

def draw_edges(corners, inner_lines=True, update_edges=False, color=PURPLE):
        corner_list = list(corners)
        sides = []

        if update_edges == True:
            edges = VGroup()
            for i in range(len(corners)):
                for j in range(i+1, len(corners)):
                    edge =  Line(
                                color = color,
                                start = corners[i].get_center(),
                                end = corners[j].get_center()
                            )
                    edge.add_updater(
                        lambda line, i=i, j=j: line.put_start_and_end_on(
                            corners[i].get_center(),
                            corners[j].get_center()
                        )
                    )
                    edges.add(edge)
            return edges

        for i in range(len(corner_list)):
            for j in range(i+1, len(corner_list)):
                sides.append(Group(corner_list[i], corner_list[j]))
        if inner_lines == False:
            sides.pop(2)
            sides.pop(2)

        edges_list = []
        for side in sides:
            edges_list.append(Line(
                stroke_width = 5,
                color = color
            ).put_start_and_end_on(side[0].get_center(), side[1].get_center())
            )


        edges = Group()
        for edge in edges_list:
            edges.add(edge)
                
        return edges

def rotation_matrix(angle, degrees=True):
    if degrees == True:
        angle = angle * np.pi/180
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0,             0,              0]
    ])
    return rotation_matrix

class TrueIntroduction(Scene):
    def construct(self):
        statement = MathTex("\\textbf{Bloch's Theorem}").scale(2)
        self.add(statement)

class Introduction(ThreeDScene):
    def construct(self):
        lattice_points_1D = VGroup()
        ah_lattice = 1 #Lattice constant horizontal
        av_lattice = 1 #Lattice constant vertical
        point_count_h = 16 #Horizontal point count

        for i in range(point_count_h):
            lattice_points_1D.add(Dot(color=BLUE).move_to(LEFT*10+RIGHT*ah_lattice*i))
        lattice_points_1D = lattice_points_1D.move_to([0,0,0])

        self.play(Write(lattice_points_1D), run_time=2)
        infinite_arrows = VGroup(
            Arrow(start=[-1,0,0], end=[-4,0,0]),
            Arrow(start=[1,0,0], end=[4,0,0])
        ).shift(UP*1)
        infinite_arrows.add(MathTex("\\infty").next_to(infinite_arrows[0].get_end(), LEFT))
        infinite_arrows.add(MathTex("\\infty").next_to(infinite_arrows[1].get_end(), RIGHT))

        self.wait() 
        self.play(Write(infinite_arrows))
        self.wait()
        self.play(lattice_points_1D.animate.shift(RIGHT*1.5*ah_lattice), run_time=3)
        self.wait()
        self.play(lattice_points_1D.animate.shift(RIGHT*1.5*ah_lattice), run_time=1.5)
        self.play(Unwrite(infinite_arrows))
        self.wait()

        point_count_v = 17 #Vertical point coint
        lattice_points_2D = VGroup()
        for i in range(point_count_v):
            lattice_points_2D.add(
                lattice_points_1D.copy().shift(UP*4.5+DOWN*av_lattice*i)
            )
        
        lattice_points_2D.move_to([0,0,0])
        self.play(Write(lattice_points_2D), FadeOut(lattice_points_1D), run_time=4)
        self.remove(lattice_points_1D)
        self.wait()

        lattice_constant_val = ValueTracker(np.linalg.norm(lattice_points_2D[0][0].get_center() - lattice_points_2D[0][1].get_center())) 
        lattice_constant = VGroup(
            MathTex("a"),
            DecimalNumber(lattice_constant_val.get_value(), num_decimal_places=3, include_sign=False, unit=None)
        ).move_to([4.2,2,0]).set_color(RED).scale(1.5)
        lattice_constant[1] = lattice_constant[1].next_to(lattice_constant[0], DOWN, buff=0.2)

        lattice_constant[1].add_updater(
            lambda d: d.set_value(lattice_constant_val.get_value())
        )

        lattice_bracket = BraceBetweenPoints(lattice_points_2D[10][10].get_center(), lattice_points_2D[10][11].get_center()).set_color(RED)
    
        self.play(Write(lattice_constant), DrawBorderThenFill(lattice_bracket))

        def track_lattice_constant(scale_factor):
            return (lattice_constant_val.animate.set_value(np.linalg.norm(lattice_points_2D.copy().scale(scale_factor)[0][0].get_center() 
                                                                          - lattice_points_2D.copy().scale(scale_factor)[0][1].get_center())))

        self.wait(2)
        self.play(FadeOut(lattice_bracket))
        self.play(
            lattice_points_2D.animate.scale(0.9), 
            track_lattice_constant(0.9),
            run_time=1.5,
        )
        self.wait()
        self.play(Rotate(lattice_points_2D, about_point=ORIGIN), run_time=5)
        self.play(Unwrite(lattice_constant))
        self.wait()
        self.move_camera(
            phi = -100*DEGREES,
            theta = 70*DEGREES,
            zoom = 0.25,
            frame_center = [0,0,-5],
            run_time=2
        )
        
        self.wait()
        threeD_lattice_layer = VGroup()
        for i in range(len(lattice_points_2D)):
            for j in range(len(lattice_points_2D[0])):
                threeD_lattice_layer.add(Dot().move_to(lattice_points_2D[i][j].get_center()).scale(1).set_color(BLUE))
                
        
        self.play(Create(threeD_lattice_layer))

        threeD_lattice = VGroup()
        for i in range(12):
            threeD_lattice.add(threeD_lattice_layer.copy().shift(np.array([0,0,-1])*i*lattice_constant_val.get_value()))

        self.wait()

        self.play(Create(threeD_lattice))

        self.wait()
        self.move_camera(
            phi = -120*DEGREES,
            run_time = 6
        )

        self.wait()

        self.move_camera(
            zoom = 4
        )

        self.wait()
        self.move_camera(
            theta = 50*DEGREES,
            phi = -100*DEGREES,
            rum_time = 7
        )
        self.wait()

        ghost_lattice = Group(*[mob.copy().set_opacity(0.2) for mob in self.mobjects])
        not_ghost = Group(*[mob for mob in self.mobjects])

        self.add(ghost_lattice)
        self.remove(not_ghost)
        self.play(
            *[mob.animate.shift(([0, 0, lattice_constant_val.get_value()])) for mob in not_ghost],
            run_time = 5
        )
        self.play(
            *[mob.animate.shift(([0 ,-lattice_constant_val.get_value(), 0])) for mob in not_ghost],
            run_time = 5
        )
        self.wait()
        #Pointing out finiteness of lattice
        cam_phi, cam_theta = self.camera.phi, self.camera.theta
        cam_frame_center = self.camera.frame_center
        relative_position = cam_frame_center -  np.array([
            np.sin(cam_phi) * np.cos(cam_theta),
            np.sin(cam_phi) * np.sin(cam_theta),
            np.cos(cam_phi)
        ]) 

        self.wait()
        self.remove(ghost_lattice)

        self.wait()

        self.move_camera(
            frame_center = [0,0,0],
            phi = 0*DEGREES,
            theta = -90*DEGREES,
            zoom = 1,
            run_time = 5
        )
        self.wait()
        self.play(
            FadeOut(threeD_lattice),
        )
        self.play(FadeIn(lattice_points_1D.move_to([0,0,0])), FadeOut(lattice_points_2D), FadeOut(threeD_lattice_layer))
        self.wait()

class PeriodicBoundaryConditions(Scene):
    def construct(self):
        lattice_points_1D = VGroup()
        ah_lattice = 1 #Lattice constant horizontal
        av_lattice = 1 #Lattice constant vertical
        point_count_h = 20 #Horizontal point count

        for i in range(point_count_h):
            lattice_points_1D.add(Dot(color=BLUE).move_to(LEFT*(10+2*ah_lattice)+RIGHT*ah_lattice*i))
        lattice_points_1D = lattice_points_1D.move_to([0,0,0])

        self.play(FadeIn(lattice_points_1D))
        self.wait()

        lattice_bracket = always_redraw(
            lambda: BraceBetweenPoints(
                lattice_points_1D[int(len(lattice_points_1D)/2)-1].get_center(),
                lattice_points_1D[int(len(lattice_points_1D)/2)].get_center()
                ).set_color(RED)
        )

        self.play(Write(lattice_bracket))
        self.wait()

        lattice_constant_val = ValueTracker(np.linalg.norm(lattice_points_1D[0].get_center() - lattice_points_1D[1].get_center())) 

        lattice_constant = VGroup(
            MathTex("a"),
            DecimalNumber(lattice_constant_val.get_value(), num_decimal_places=3, include_sign=False, unit=None)
        ).move_to([0,.8,0]).set_color(RED).scale(1.5)
        lattice_constant[1] = lattice_constant[1].next_to(lattice_constant[0], DOWN, buff=0.2)        
        lattice_constant[1].add_updater(
            lambda d: d.set_value(lattice_constant_val.get_value())
        )
        lattice_constant.add(SurroundingRectangle(lattice_constant, color=DARK_BLUE))
        
        self.play(lattice_points_1D.animate.shift(DOWN*1))
        self.wait()
        self.play(Write(lattice_constant))
        self.wait()
        self.play(lattice_points_1D.animate.scale(1.5), lattice_constant_val.animate.set_value(1.5))
        self.wait()
        self.play(lattice_points_1D.animate.scale(0.8), lattice_constant_val.animate.set_value(1.5*0.8))
        self.wait()
        new_lb = always_redraw(lambda: lattice_bracket.copy().shift(RIGHT*lattice_constant_val.get_value()))
        self.play(Write(new_lb))
        self.wait()
        self.play(lattice_points_1D.animate.shift(LEFT*lattice_constant_val.get_value()))
        self.remove(lattice_bracket)
        self.wait()
        self.play(
            lattice_constant.animate.next_to(new_lb.get_tip(), DOWN),
        )

        #Adding sine wave
        wavelength = 2
        ax = Axes(
            x_range = [-len(lattice_points_1D)/2+1, len(lattice_points_1D)/2],
            y_range = [-2,2],
            x_length = lattice_points_1D[-1].get_center()[0] - lattice_points_1D[0].get_center()[0],
            y_length = 4
        ).shift(DOWN*1)
        wave_func = ax.plot(lambda x: np.cos(2*PI/wavelength * x)**2, color=BLUE_E)

        self.play(Create(wave_func), run_time=2)

        self.wait()
        self.play(wave_func.animate.shift(RIGHT*lattice_constant_val.get_value()), run_time=2)
        self.wait()

        BVK = MathTex("\\textbf{Translational Symmetry}").shift(UP*1.5).scale(1.2)

        self.play(Write(BVK))
        self.wait()
        self.play(
            BVK.animate.scale(0.8).move_to([0,3,0]),
        )

        PSI = MathTex(
            "|",                                #0
            "\\psi(x)",                            #1
            "|^2",                              #2
            "=",                       #3     
            "|",                         #4
            "\\psi(x+sa)",                #5
            "|^2",                     #6
            ", s = 0, 1, 2... ",        #7
            # "\\text{ Number of atoms in crystal}" #8
        ).move_to([0,1.5,0]).scale(0.9)
        
        self.play(Write(PSI[1]))
        self.wait()
        self.play(Write(PSI[3]), Write(PSI[5]))
        self.wait()
        self.play(Write(PSI[0]), Write(PSI[2]), Write(PSI[4]), Write(PSI[6]))
        self.wait()
        self.play(Write(PSI[7]))
        self.wait() 
            
        hopping_dot = Dot().move_to(lattice_points_1D[8].get_center()+[0,1,0]).set_color(GREEN)
        self.play(FadeIn(hopping_dot))
        self.play(hopping_dot.animate.shift(RIGHT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(RIGHT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(2*RIGHT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(LEFT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(2*LEFT*lattice_constant_val.get_value()))
        self.play(FadeOut(hopping_dot))
        self.wait()       

        # self.play(FadeIn(ax.shift(LEFT*3*lattice_constant_val.get_value())))
        complex_exponential = MathTex(
            "|C| = 1",   #0
            "\\implies", #1
            "C",         #2
            "=",         #3
            "e^",        #4
            "{i",        #5
            "\\theta}"         #6
        )

        nth_roots = MathTex(
            "C",      #0
            "C",          #1
            "= e^{2\\pi is/N}",     #2
            ", s = 0, 1, 2... N-1"    #3
        )
        temp = MathTex("C").move_to(PSI[5].get_center()+LEFT*0.3)
        self.wait()
        PSI_GROUP = VGroup(PSI, nth_roots)
        self.play(
            FadeIn(nth_roots[0].move_to(PSI[1].get_center()+LEFT*0.2)),
            *[mob.animate.shift(RIGHT*0.5) for mob in PSI[1:]]
        )
        self.wait()
        self.play(FadeIn(complex_exponential[0].shift(UP*0.5)))
        """
        self.play(
            FadeIn(temp),
            *[mob.animate.shift(RIGHT*0.5) for mob in PSI[5:]]
        )
        self.wait(2)
        self.play(
            *[mob.animate.shift(LEFT*0.5) for mob in PSI[5:]],
            FadeOut(temp)
        )
        """
        self.wait(2)

        self.play(
            *[
                Unwrite(PSI[i]) for i in np.arange(0, 8, 2)
            ]
        )
        self.play(
            PSI[1].animate.shift(RIGHT*0.4),
            nth_roots[0].animate.shift(RIGHT*0.4),
            PSI[3].animate.shift(RIGHT*0.1)
        )

        self.wait(2)

        self.play(FadeIn(complex_exponential[1:].shift(UP*0.5)))
        self.wait()
        self.play(Circumscribe(complex_exponential[6]))
        self.wait()

        r = 2.5 #circle radius
        ghost_circle = Circle(radius=r).set_opacity(0)
        circular_crystal = VGroup()
        N = 12 #Number of points on circle
        for i in range(N):
            angle = 2*PI/N * i
            circular_crystal.add(
                Dot(color=BLUE).move_to(ghost_circle.get_center() + r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        temp = self.mobjects
        self.play(
            ReplacementTransform(lattice_points_1D.copy(), circular_crystal),
            *[FadeOut(mob) for mob in self.mobjects]
        )
        self.wait()
        self.play(*[FadeIn(mob) for mob in temp], FadeOut(circular_crystal))
        self.wait()

        self.play(
            wave_func.animate.shift(DOWN*2.2),
            lattice_points_1D.animate.shift(DOWN*2.2),
            Unwrite(lattice_constant)
        )
        self.wait()
        self.play(
            ReplacementTransform(complex_exponential, nth_roots[1:3].shift(UP*0.5))
        )
        self.wait()
        self.play(
            Write(nth_roots[3].shift(UP*0.5))
        )
        self.wait(2)
        bloch_theorem = MathTex("\\psi(x+sa) = e^{","i2\\pi s/N","}\\psi(x)").move_to([0,-0.5,0])
        self.play(Write(bloch_theorem))
        self.wait(2)
        bloch_title = Tex("Bloch's Theorem").next_to(bloch_theorem, DOWN, buff=0.3)
        bloch_theorem = VGroup(bloch_theorem, bloch_title)
        bloch_rectangle = SurroundingRectangle(bloch_theorem, color=GOLD)
        bloch_theorem.add(bloch_rectangle)
        self.play(FadeIn(bloch_title), FadeIn(bloch_rectangle))
        self.wait()

        #Adding 'k'
        k = MathTex("^{ika}")
        k_expression = MathTex("k =", "2\\pi s/Na")
        bloch_theorem.add(k, k_expression)
        self.play(
            ReplacementTransform(bloch_theorem[0][1], k.move_to(bloch_theorem[0][1].get_center()+LEFT*0.3)),
            bloch_theorem[0][2].animate.shift(LEFT*0.6)
        )
        self.wait()
        new_rectangle = SurroundingRectangle(bloch_theorem, color=GOLD).move_to(bloch_title.get_center()+LEFT*0.2+UP*0.2)
        self.play(
            Write(k_expression.move_to(bloch_title.get_center()+UP*0.2)), 
            bloch_title.animate.shift(DOWN*0.3),
            ReplacementTransform(bloch_rectangle, new_rectangle)
        )
        self.wait()
        bloch_1D = Tex("Bloch's Theorem", " (In 1D)").move_to(bloch_title.get_center()+LEFT*0.2)
        self.play(
            bloch_title.animate.move_to(bloch_1D[0].get_center()),
            FadeIn(bloch_1D[1])
        )
        bloch_theorem.add(bloch_1D)

        self.wait(3)
        no_fade_group = VGroup(
            bloch_theorem, wave_func, lattice_points_1D, lattice_bracket
        )
        fade_group = Group()
        for mob in self.mobjects:
            if mob in no_fade_group or mob in bloch_theorem or mob is new_lb or mob in bloch_1D or mob is new_rectangle:
                continue
            fade_group.add(mob)
        self.play(
            *[
                Unwrite(mob) for mob in VGroup(nth_roots, complex_exponential, PSI[1],
                                               PSI[3], PSI[5], PSI[7], BVK)
            ],
            run_time = 0.7
        )
        self.wait()
        #Final form of bloch's theorem in 1D
        satisfactoy_wf = MathTex(
            "\\psi(x)",     #0
            "=",            #1
            "e^{ikx}",      #2
            "\\implies",    #3
            "\\psi(x+",     #4
            "s",            #5
            "a)",           #6
            "=",            #7
            "e^{i2\\pi ",   #8
            "s",            #9
            "/N}",          #10             
            "e^{ikx}",      #11   
        ).shift(UP*2)
        psi_replacement = MathTex("N", "^N", "\\psi(x)")
        self.play(Write(satisfactoy_wf[0:3]))
        self.wait()
        self.play(Write(satisfactoy_wf[3:]))
        self.wait()
        self.play(
            ReplacementTransform(satisfactoy_wf[5], psi_replacement[0].move_to(satisfactoy_wf[5].get_center()+UP*0.05+LEFT*0.03)),
            ReplacementTransform(satisfactoy_wf[9], psi_replacement[1].move_to(satisfactoy_wf[9].get_center()+UP*0.03))
        )
        self.wait()
        self.play(
            ReplacementTransform(satisfactoy_wf[11], psi_replacement[2].move_to(satisfactoy_wf[11].get_center()+RIGHT*0.15))
        )
        self.wait()
        self.play(
            FadeOut(satisfactoy_wf[8]), FadeOut(satisfactoy_wf[10]), FadeOut(psi_replacement[1]),
            psi_replacement[2].animate.shift(LEFT*1.55+DOWN*0.05)
        )
        self.wait()

        true_psi = MathTex(
            "\\psi(x)",     #0
            "=",            #1
            "e^{ikx}",      #2
            "u",            #3
            "(x)",          #4
            ",",            #5
            "u(x+a)",       #6
            "=",            #7
            "u(x)"          #8
        ).shift(UP*0.8)
        self.play(Write(true_psi[0:5]))
        self.wait()
        self.play(
            Write(true_psi[5]),
            Write(true_psi[6:].shift(RIGHT*0.2))
        )
        self.wait()
        self.play(
            Circumscribe(true_psi[3:]),
            run_time = 2
        )
        self.play(
            Circumscribe(Rectangle(width=13, height=1).move_to([0, wave_func.get_center()[1], 0])),
            run_time = 2
        )
        self.wait()
        final_rect = SurroundingRectangle(
            VGroup(new_rectangle.copy(), true_psi)
        ).set_color(GOLD).shift(UP*0.15)
        self.play(
            FadeOut(bloch_rectangle),
            ReplacementTransform(new_rectangle, final_rect)
        )
        self.play(
            k_expression.animate.shift(UP*0.25),
            bloch_theorem[0][0].animate.shift(UP*0.5),
            bloch_theorem[0][2].animate.shift(UP*0.5),
            k.animate.shift(UP*0.5)
        )

        BLOCH_THEOREM = VGroup(
            k_expression.copy(), bloch_theorem[0][0].copy(), bloch_theorem[0][2].copy(), k.copy(),
            final_rect.copy(), bloch_title.copy(), bloch_1D[1].copy(), true_psi.copy()
        )
        self.wait()
        self.play(
            *[
                Unwrite(mob) for mob in VGroup(
                    k_expression, bloch_theorem[0][0], bloch_theorem[0][2], k, final_rect, psi_replacement[0], psi_replacement[2],
                    bloch_title, bloch_1D[1], true_psi, wave_func, lattice_points_1D, new_lb, satisfactoy_wf[0:5], satisfactoy_wf[6:8]
                )
            ],
            run_time = 0.7
        )
        self.wait()
        self.play(
            FadeIn(BLOCH_THEOREM.shift(UP*2))
        )
        self.wait()

        takeaways = BulletedList(
            "Periodic unit cell function/potential",
            "Shifting by lattice constant picks up a phase factor",
            "Crystal momentum 'k'; reciprocal space"
        ).shift(DOWN*1.8)

        self.play(Write(takeaways[0]))
        self.wait()
        
        self.play(Write(takeaways[1]))
        self.wait()

        self.play(
            Circumscribe(BLOCH_THEOREM[0])
        )
        self.wait()
        
        standard_k = MathTex("k = 2\\pi / \\lambda").move_to([0,takeaways[2].get_center()[1],0])
        self.play(Write(standard_k))
        self.wait()
        self.play(ReplacementTransform(standard_k, takeaways[2]))
        self.wait()
        
        
        """
        #Begin Schrodinger
        periodic_potential = MathTex("U(x + a) = U(x)").move_to([-3,2,0]).scale(0.9)
        schrodinger_eqn = MathTex(
            "\\hat{H}\\psi = E\\psi",                                                                      #0
            "\\left[-\\frac{\\hbar^2}{2 m}\\frac{\\partial^2}{\\partial x^2} + U(x)\\right]\\psi = E\\psi" #1
        ).move_to([4,2,0]).scale(0.9)
    
        self.play(Write(schrodinger_eqn[0]))
        self.wait()
        self.play(Write(schrodinger_eqn[1].shift(DOWN*1+LEFT*3)))
        self.wait()
        self.play(Write(periodic_potential))
        self.wait()
        """
        self.wait()
        
class PhaseFactor(ThreeDScene):
    def construct(self):
        # Define the parametric equations for a spiral
        self.set_camera_orientation(phi=0*DEGREES, theta=-90*DEGREES)

        axes = ThreeDAxes(
        ).scale(2)
        # axes_labels = axes.get_axis_labels(x_label="Im", y_label="x", z_label = 'z')

        def unit_cell_func(t):
            x = np.cos(6*t)**2
            y = 0 #np.sin(2*t)**2
            z = t
            return axes.c2p(z, x, y)
        
        def phase_func(t):
            x = np.cos(t)
            y = np.sin(t)
            z = t
            return axes.c2p(z, x, y)
        
        def wave_func(t):
            x = np.cos(6*t)**2 * np.cos(t)
            y = np.sin(t) * np.cos(6*t)**2
            z = t
            return axes.c2p(z, x, y)
        
        def wave_function(t):
            x = t
            y = np.sin(t) * np.cos(6*t)**2
            z = np.cos(6*t)**2 * np.cos(t)
            return axes.c2p(x, y, z)
        
        unit_cell = ParametricFunction(
            unit_cell_func,
            t_range = [-np.pi, np.pi],
            color = PURPLE
        )

        phase_factor = ParametricFunction(
            phase_func,
            t_range = [-np.pi, np.pi],
            color = RED
        )

        wave_curve = ParametricFunction(
            wave_func,
            t_range = [-np.pi, np.pi],
            color = GOLD
        )

        axes_labels = VGroup(
            always_redraw(lambda: MathTex("x").move_to(axes.c2p(3.8, 0.2, 0)).scale(0.7)),
            always_redraw(lambda: Tex("Re").move_to(axes.c2p(0.2, 1.7, 0)).scale(0.7)),
            always_redraw(lambda: Tex("Im").move_to(axes.c2p(0, 0.2, 3.4)).scale(0.7)),
        )

        self.play(Create(axes), Create(axes_labels))
        uc = always_redraw(lambda: MathTex("u(x)").move_to(axes.c2p(-3,1.2,0)).set_color(PURPLE))
        phase_eqn = always_redraw(lambda: MathTex("e^{ikx}").move_to(axes.c2p(-0.6,1.2,-0.4)).set_color(RED))
        wave_eqn = always_redraw(lambda: MathTex("\\psi (x)").move_to(axes.c2p(0.6,1.2,0.4)).set_color(GOLD))

        self.play(Create(unit_cell), FadeIn(uc))
        self.wait()
        self.play(Create(phase_factor), FadeIn(phase_eqn))
        self.wait()

        self.move_camera(
            phi = -120*DEGREES,
            run_time = 5
        )
        self.wait()
        self.move_camera(
            phi = -10*DEGREES,
            run_time = 5
        )
        self.wait()
        self.play(ReplacementTransform(unit_cell, wave_curve), FadeIn(wave_eqn))
        self.wait()
        new_uc = ParametricFunction(
            unit_cell_func,
            t_range = [-np.pi, np.pi],
            color = PURPLE
        )
        self.wait()
        self.move_camera(
            phi = -120*DEGREES,
            run_time = 5
        )
        self.wait()
        self.move_camera(
            phi = -10*DEGREES,
            run_time = 5
        )
        self.wait()

        self.play(
            *[
                Rotate(mob, axis = UP, angle = -60*DEGREES) for mob in VGroup(axes, wave_curve, phase_factor)
            ],
            run_time = 12
        )
        self.wait()
        uc_cut = ParametricFunction(
            unit_cell_func,
            t_range = [np.pi/3, 5*np.pi/6],
            color = PURPLE
        )

        wave_curve_cut = ParametricFunction(
            wave_func,
            t_range = [np.pi/3, 5*np.pi/6],
            color = PURPLE
        )
        self.play(FadeIn(uc_cut))
        self.play(
            ReplacementTransform(uc_cut, wave_curve_cut)
        )
        self.wait(2)
        self.play(FadeOut(VGroup(wave_curve_cut, uc_cut)))
        self.play(
            *[
                Rotate(mob, axis = UP, angle = 60*DEGREES) for mob in VGroup(axes, wave_curve, phase_factor)
            ],
            run_time = 12
        )
        self.wait()
        traveller = Dot(color = BLUE_D).move_to(wave_func(0)).scale(2)
        distance_lines = VGroup()
        time = ValueTracker(0)
        def traveller_updater(axes, vt):
            def updater(mob):
                mob.move_to(axes.c2p(
                    vt.get_value(),
                    np.cos(6*vt.get_value())**2 * np.cos(vt.get_value()),
                    np.cos(6*vt.get_value())**2 * np.sin(vt.get_value())
                ))
            return updater
        traveller.add_updater(traveller_updater(axes, time))
        self.play(DrawBorderThenFill(traveller), run_time = 2)
        self.wait()
       
        def add_distance_lines(group, func, vals, vt, color = BLUE_D, run_time=1):
            lines = []
            for i, value in enumerate(vals):
                self.play(vt.animate.set_value(vals[i]), run_time=run_time)
                lines.append(always_redraw(lambda i=i: Line(
                        start = func(vals[i]),
                        end = axes.c2p(vals[i], 0, 0),
                        color = color
                    )))
                
                group.add(lines[i])
                self.play(Write(lines[i]), run_time=run_time)

        self.play(time.animate.set_value(-1)) 
        distance_lines.add(
            always_redraw(lambda: Line(
                start = wave_func(-1),
                end = axes.c2p(-1, 0, 0),
                color = BLUE_D
            ))
        )
        self.play(Write(distance_lines[0]))
        
        self.play(
            *[
                Rotate(mob, axis = UP, angle = 40*DEGREES) for mob in VGroup(axes, wave_curve, phase_factor)
            ],
            run_time = 2
        )

        add_distance_lines(
            distance_lines,
            wave_func,
            [np.pi/6-1, np.pi/3-1, np.pi/2-1, 2*np.pi/3-1, 5*np.pi/6-1],
            time,
            run_time = 0.5
        )
        
        self.wait()
        self.play(
            *[
                mob.animate.rotate(axis=UP, angle=-40*DEGREES).shift(UP*1) for mob in VGroup(axes, wave_curve, phase_factor)
            ],
        )
        self.play(Unwrite(axes), Unwrite(axes_labels))
        self.wait()
        new_distance_lines = VGroup()
        for mob in distance_lines:
            new_distance_lines.add(mob.copy())
        self.add(new_distance_lines)

        self.play(
            *[
                new_distance_lines[i].animate.put_start_and_end_on([-1, -i*0.5, 0], [1, -i*0.5, 0]) for i in range(len(new_distance_lines))
            ]
        )
        self.play(Unwrite(distance_lines), run_time = 2)
        self.wait()
        self.play(Unwrite(new_distance_lines), run_time = 2)
        self.wait()

        unallowed_points = VGroup()
        for i in range(12):
            unallowed_points.add(
                Dot().move_to(
                    axes.c2p(-11*np.pi/12+i*np.pi/6, 0, 0)
                )
            )
        self.play(Write(unallowed_points))
        self.wait()
        Unallowed_Energies = Tex("Unallowed Energies").shift(DOWN*0.5).scale(1.2)
        self.play(Write(Unallowed_Energies), run_time=1)
        self.wait()
        self.play(
            *[
                Unwrite(mob) for mob in VGroup(unallowed_points, wave_curve, phase_factor, traveller, uc, phase_eqn, wave_eqn, Unallowed_Energies)
            ],
            run_time = 1
        )

        self.wait()
        
        
class BornvonKarman(Scene):
    def construct(self):

        r = 2.5 #circle radius
        ghost_circle = Circle(radius=r).set_opacity(0)
        circular_crystal = VGroup()
        N = 12 #Number of points on circle
        for i in range(N):
            angle = 2*PI/N * i
            circular_crystal.add(
                Dot(color=BLUE).move_to(ghost_circle.get_center() + r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        return_circular_crystal = circular_crystal.copy() #What to finish this scene with
        self.add(circular_crystal)
        self.wait()

        traveller = Dot(color=RED).scale(1.3)
        trail = TracedPath(traveller.get_center, dissipating_time=0.3, stroke_opacity=[0, 1], stroke_width=8).set_color(RED)
        self.add(traveller, trail)

        for point in circular_crystal:
            angle = 2*PI/N
            self.play(
                traveller.animate(path_arc=angle).move_to(point.get_center()),
                run_time = 0.65
        )
        self.play(
            traveller.animate(path_arc=2*PI/N).move_to(circular_crystal[0].get_center()),
            run_time = 0.65
        )
        self.wait()
        self.remove(trail)
        BVK = Text("Periodic Boundary \n      Conditions").scale(0.8).move_to([0,-2.3, 0])
        self.play(
            circular_crystal.animate.move_to(ORIGIN+r*np.array([0,1,0])),
            traveller.animate.move_to(circular_crystal[N-3].get_center()+r*np.array([0,1,0])),
            Write(BVK)
        )
        crystal_trail = TracedPath(circular_crystal[N-3].get_center, dissipating_time=30, stroke_opacity=[0,2], stroke_width=8).set_color(RED)
        self.add(crystal_trail)
        for i in range(N):
            angle = 2*PI/N
            self.play(
                traveller.animate.move_to(circular_crystal[N-i-3].get_center()),
                Rotate(circular_crystal, angle=angle),
                run_time=0.65
            )
        self.wait()
        self.play(
            FadeOut(traveller), Unwrite(crystal_trail)
        )
        self.wait()
    
        r_tracker = ValueTracker(2.5)
        radius_tracker = DecimalNumber(r_tracker.get_value(),num_decimal_places=2, include_sign=False, unit=None).shift(DOWN*1)
        radius = Tex("Radius:").next_to(radius_tracker, LEFT, buff=0.2)
        radius_tracking_group = VGroup(radius_tracker, radius)

        radius_tracker.add_updater(
            lambda d: d.set_value(r_tracker.get_value())
        )
        self.play(
            FadeIn(radius_tracking_group)
        )
        self.wait()


        big_r = 20
        big_N = 250
        big_circle = VGroup()
        for i in range(big_N):
            angle = 2*PI/big_N * i
            big_circle.add(
                Dot(color=BLUE).move_to([0,big_r,0] + big_r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        self.play(
            r_tracker.animate.set_value(20),
            ReplacementTransform(circular_crystal, big_circle)
        )
        self.wait()

        super_big_r = 300
        super_big_N = 3000
        super_big_circle = VGroup()
        for i in range(super_big_N):
            angle = 2*PI/super_big_N * i
            super_big_circle.add(
                Dot(color=BLUE).move_to([0,super_big_r,0] + super_big_r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        self.play(
            r_tracker.animate.set_value(300),
            ReplacementTransform(big_circle, super_big_circle)
        )
        self.wait()
        straight_line = Line(start=LEFT*10, end=RIGHT*10).set_opacity(0.5)
        self.play(FadeIn(straight_line), run_time=2)
        self.wait()

        ultra_big_r = 1000
        ultra_big_N = 7000
        ultra_big_circle = VGroup()
        for i in range(ultra_big_N):
            angle = 2*PI/ultra_big_N * i
            ultra_big_circle.add(
                Dot(color=BLUE).move_to([0,ultra_big_r,0] + ultra_big_r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        self.play(
            r_tracker.animate.set_value(1000),
            ReplacementTransform(super_big_circle, ultra_big_circle),
            run_time = 6,
        )
        infty = MathTex("\\infty").move_to(radius_tracker).scale(1.5)
        self.play(Unwrite(radius_tracker), FadeIn(infty.shift(LEFT*0.3)))
        radius_tracking_group.remove(radius_tracker)
        radius_tracking_group.add(infty)

        self.wait(2)
        self.play(Unwrite(straight_line), Unwrite(radius_tracking_group))
        self.wait()
        self.play(Unwrite(BVK))
        self.play(ReplacementTransform(ultra_big_circle, return_circular_crystal), run_time=4)
        self.wait()

class NthRootofUnity(Scene):
    def construct(self):
        r = 2.5 #circle radius
        ghost_circle = Circle(radius=r).set_opacity(0)
        circular_crystal = VGroup()
        N = 12 #Number of points on circle
        for i in range(N):
            angle = 2*PI/N * i
            circular_crystal.add(
                Dot(color=BLUE).move_to(ghost_circle.get_center() + r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        return_circular_crystal = circular_crystal.copy() #What to finish this scene with
        self.add(circular_crystal)
        self.wait()
        plane = NumberPlane(
            x_range = [-3, 3],
            y_range =  [-3, 3],
            x_length = 6*r,
            y_length = 6*r
        ).set_opacity(0.3).set_color(GOLD)
        ax = Axes(
            x_range = [-2.8, 2.8],
            y_range =  [-1.6, 1.6],
            x_length = 5.6*r,
            y_length = 3.2*r
        ).set_opacity(0.4).set_color(WHITE)

        ax_labels = VGroup(
            Tex("Re").move_to([6.7,0.5,0]),
            Tex("Im").move_to([-0.8,3.7,0])
        )

        axes = VGroup(ax, plane, *[mob for mob in ax_labels])
        self.play(FadeIn(axes))
        self.wait(2)

        traveller = Dot(color=RED).scale(1.3)
        trail = TracedPath(traveller.get_center, dissipating_time=0.3, stroke_opacity=[0, 1], stroke_width=8).set_color(RED)
        self.add(traveller, trail)

        nth_roots = VGroup(
            MathTex("e^{i0\\pi/6"),
            MathTex("e^{i\\pi/6}"),
            MathTex("e^{i\\pi/3}"),
            MathTex("e^{i\\pi/2}"),
            MathTex("e^{i2\\pi/3}"),
            MathTex("e^{i5\\pi/6}"),
            MathTex("e^{i\\pi}"),
            MathTex("e^{i7\\pi/6}"),
            MathTex("e^{i4\\pi/3}"),
            MathTex("e^{i3\\pi/2}"),
            MathTex("e^{i5\\pi/3}"),
            MathTex("e^{i11\\pi/6}"),
        )


        angle = 2*PI/N
        for i, point in enumerate(circular_crystal):
            angle = 2*PI/N
            self.play(
                traveller.animate(path_arc=angle).move_to(point.get_center()),
                FadeIn(nth_roots[i].move_to(point.get_center()*1.2+RIGHT*0.2)),
                run_time=0.8
            )
        self.play(
            traveller.animate(path_arc=angle).move_to(circular_crystal[0].get_center())
        )
        self.wait()
        self.remove(trail)
        self.wait()
        nth_root_forms = VGroup(
            MathTex("e^{i2\\pi(0)/12}"),
            MathTex("e^{i2\\pi(1)/12}"),
            MathTex("e^{i2\\pi(2)/12}"),
            MathTex("e^{i2\\pi(3)/12}"),
            MathTex("e^{i2\\pi(4)/12}"),
            MathTex("e^{i2\\pi(5)/12}"),
            MathTex("e^{i2\\pi(6)/12}"),
            MathTex("e^{i2\\pi(7)/12}"),
            MathTex("e^{i2\\pi(8)/12}"),
            MathTex("e^{i2\\pi(9)/12}"),
            MathTex("e^{i2\\pi(10)/12}"),
            MathTex("e^{i2\\pi(11)/12}"),
        )
        for i, point in enumerate(circular_crystal):
            angle = 2*PI/N
            nth_root_forms[i].scale(0.8).move_to(point.get_center()*1.2+RIGHT*0.1)

        self.play(
            *[ReplacementTransform(nth_roots[i], nth_root_forms[i]) for i, root in enumerate(nth_roots)]
        )
        self.wait(2)
        nth_root_statement = Tex("""Nth Roots\n
                                 Of Unity""")
        nth_root_box = BackgroundRectangle(nth_root_statement)
        nth_root_rect = SurroundingRectangle(nth_root_statement).set_color(GOLD)
        nth_root_statement_group = VGroup(nth_root_box, nth_root_statement, nth_root_rect)
        self.play(FadeIn(nth_root_statement_group))
        self.wait()

        nth_root_form = MathTex("e^{i2\\pi s/N}").shift(DOWN*0.55)
        nth_root_form_box = BackgroundRectangle(nth_root_form).set_opacity(0.8)
        nth_root_form_group = VGroup(nth_root_form_box, nth_root_form)
        self.play(
            FadeIn(nth_root_form_group),
            nth_root_statement_group.animate.shift(UP*0.3),
            ReplacementTransform(nth_root_rect, SurroundingRectangle(VGroup(nth_root_statement.copy().shift(UP*.3), nth_root_form)).set_color(GOLD)),
            ReplacementTransform(nth_root_box, BackgroundRectangle(VGroup(nth_root_statement.copy().shift(UP*.3), nth_root_form)).set_color(BLACK))
        )
        self.wait()


#ADDITIONS FOR VIDEO (CLIPS THAT I THOUGHT NECESSARY, BUT HAD ALREADY RENDERED THE SCENE)
class Long1D(Scene):
    """Stretches short line of dots into long line,
    then, comes back
    """
    def construct(self):
        lattice_points_1D = VGroup()
        ah_lattice = 1 #Lattice constant horizontal
        point_count_h = 16 #Horizontal point count

        for i in range(point_count_h):
            lattice_points_1D.add(Dot(color=BLUE).move_to(LEFT*10+RIGHT*ah_lattice*i))
        lattice_points_1D = lattice_points_1D.move_to([0,0,0])
        self.add(lattice_points_1D)
        backup = lattice_points_1D.copy()

        long_lattice = VGroup()
        new_la = 0.15
        many_points = 100
        for i in range(many_points):
            long_lattice.add(Dot(color=BLUE).move_to(LEFT*10+RIGHT*new_la*i))
        long_lattice = long_lattice.move_to([0,0,0])

        self.play(
            ReplacementTransform(lattice_points_1D, long_lattice),
            run_time = 2
        )

        self.wait()
        self.play(
            ReplacementTransform(long_lattice, backup)
        )
        self.wait()

class WaveFunctionPeriodicity(Scene):
    def construct(self):
        r = 2.5 #circle radius
        ghost_circle = Circle(radius=r).set_opacity(0)
        circular_crystal = VGroup()
        N = 12 #Number of points on circle
        for i in range(N):
            angle = 2*PI/N * i
            circular_crystal.add(
                Dot(color=BLUE).move_to(ghost_circle.get_center() + r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        return_circular_crystal = circular_crystal.copy() #What to finish this scene with
        self.add(circular_crystal)

        periodicity = MathTex("\\psi(x) = \\psi(x + Na)").shift(UP*0.5)
        nth_root = MathTex("C^N = 1").shift(DOWN*0.5)
        self.play(Write(periodicity))
        self.wait(2)
        self.play(Write(nth_root))
        self.wait()
        self.play(FadeOut(periodicity), FadeOut(nth_root))
        self.wait()

class Prerequisites(Scene):
    def construct(self):

        prerequisites = VGroup(
            MathTex("\\textbf{Useful Knowledge for This Video}").scale(1.2).shift(UP*0.75),
            BulletedList(
            "Basic Quantum Mechanics",
            "Basic Complex Analysis"
            ).scale(1).shift(DOWN*0.75)
        )
        self.add(prerequisites)

class EquationCorrection(Scene):
    def construct(self):
        lattice_points_1D = VGroup()
        ah_lattice = 1 #Lattice constant horizontal
        point_count_h = 20 #Horizontal point count

        for i in range(point_count_h):
            lattice_points_1D.add(Dot(color=BLUE).move_to(LEFT*(10+2*ah_lattice)+RIGHT*ah_lattice*i))
        lattice_points_1D = lattice_points_1D.move_to([0,0,0])

        self.play(FadeIn(lattice_points_1D))
        self.wait()

        lattice_bracket = always_redraw(
            lambda: BraceBetweenPoints(
                lattice_points_1D[int(len(lattice_points_1D)/2)-1].get_center(),
                lattice_points_1D[int(len(lattice_points_1D)/2)].get_center()
                ).set_color(RED)
        )

        self.play(Write(lattice_bracket))
        self.wait()

        lattice_constant_val = ValueTracker(np.linalg.norm(lattice_points_1D[0].get_center() - lattice_points_1D[1].get_center())) 

        lattice_constant = VGroup(
            MathTex("a"),
            DecimalNumber(lattice_constant_val.get_value(), num_decimal_places=3, include_sign=False, unit=None)
        ).move_to([0,.8,0]).set_color(RED).scale(1.5)
        lattice_constant[1] = lattice_constant[1].next_to(lattice_constant[0], DOWN, buff=0.2)        
        lattice_constant[1].add_updater(
            lambda d: d.set_value(lattice_constant_val.get_value())
        )
        lattice_constant.add(SurroundingRectangle(lattice_constant, color=DARK_BLUE))
        
        self.play(lattice_points_1D.animate.shift(DOWN*1))
        self.wait()
        self.play(Write(lattice_constant))
        self.wait()
        self.play(lattice_points_1D.animate.scale(1.5), lattice_constant_val.animate.set_value(1.5))
        self.wait()
        self.play(lattice_points_1D.animate.scale(0.8), lattice_constant_val.animate.set_value(1.5*0.8))
        self.wait()
        new_lb = always_redraw(lambda: lattice_bracket.copy().shift(RIGHT*lattice_constant_val.get_value()))
        self.play(Write(new_lb))
        self.wait()
        self.play(lattice_points_1D.animate.shift(LEFT*lattice_constant_val.get_value()))
        self.remove(lattice_bracket)
        self.wait()
        self.play(
            lattice_constant.animate.next_to(new_lb.get_tip(), DOWN),
        )

        #Adding sine wave
        wavelength = 2
        ax = Axes(
            x_range = [-len(lattice_points_1D)/2+1, len(lattice_points_1D)/2],
            y_range = [-2,2],
            x_length = lattice_points_1D[-1].get_center()[0] - lattice_points_1D[0].get_center()[0],
            y_length = 4
        ).shift(DOWN*1)
        wave_func = ax.plot(lambda x: np.cos(2*PI/wavelength * x)**2, color=BLUE_E)

        self.play(Create(wave_func), run_time=2)

        self.wait()
        self.play(wave_func.animate.shift(RIGHT*lattice_constant_val.get_value()), run_time=2)
        self.wait()

        BVK = MathTex("\\textbf{Translational Symmetry}").shift(UP*1.5).scale(1.2)

        self.play(Write(BVK))
        self.wait()
        self.play(
            BVK.animate.scale(0.8).move_to([0,3,0]),
        )

        PSI = MathTex(
            "|",                                #0
            "\\psi(x)",                            #1
            "|^2",                              #2
            "=",                       #3     
            "|",                         #4
            "\\psi(x+sa)",                #5
            "|^2",                     #6
            ", s = 0, 1, 2... ",        #7
            # "\\text{ Number of atoms in crystal}" #8
        ).move_to([0,1.5,0]).scale(0.9)
        
        self.play(Write(PSI[1]))
        self.wait()
        self.play(Write(PSI[3]), Write(PSI[5]))
        self.wait()
        self.play(Write(PSI[0]), Write(PSI[2]), Write(PSI[4]), Write(PSI[6]))
        self.wait()
        self.play(Write(PSI[7]))
        self.wait() 
            
        hopping_dot = Dot().move_to(lattice_points_1D[8].get_center()+[0,1,0]).set_color(GREEN)
        self.play(FadeIn(hopping_dot))
        self.play(hopping_dot.animate.shift(RIGHT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(RIGHT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(2*RIGHT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(LEFT*lattice_constant_val.get_value()))
        self.wait()
        self.play(hopping_dot.animate.shift(2*LEFT*lattice_constant_val.get_value()))
        self.play(FadeOut(hopping_dot))
        self.wait()       

        # self.play(FadeIn(ax.shift(LEFT*3*lattice_constant_val.get_value())))
        complex_exponential = MathTex(
            "|C| = 1",   #0
            "\\implies", #1
            "C",         #2
            "=",         #3
            "e^",        #4
            "{i",        #5
            "\\theta}"         #6
        )

        nth_roots = MathTex(
            "C",      #0
            "C",          #1
            "= e^{2\\pi is/N}",     #2
            ", s = 0, 1, 2... N-1"    #3
        )
        temp = MathTex("C").move_to(PSI[5].get_center()+LEFT*0.3)
        self.wait()
        PSI_GROUP = VGroup(PSI, nth_roots)
        self.play(
            FadeIn(nth_roots[0].move_to(PSI[1].get_center()+LEFT*0.2)),
            *[mob.animate.shift(RIGHT*0.5) for mob in PSI[1:]]
        )
        self.wait()
        self.play(FadeIn(complex_exponential[0].shift(UP*0.5)))
        """
        self.play(
            FadeIn(temp),
            *[mob.animate.shift(RIGHT*0.5) for mob in PSI[5:]]
        )
        self.wait(2)
        self.play(
            *[mob.animate.shift(LEFT*0.5) for mob in PSI[5:]],
            FadeOut(temp)
        )
        """
        self.wait(2)

        self.play(
            *[
                Unwrite(PSI[i]) for i in np.arange(0, 8, 2)
            ]
        )
        self.play(
            PSI[1].animate.shift(RIGHT*0.4),
            nth_roots[0].animate.shift(RIGHT*0.4),
            PSI[3].animate.shift(RIGHT*0.1)
        )

        self.wait(2)

        self.play(FadeIn(complex_exponential[1:].shift(UP*0.5)))
        self.wait()
        self.play(Circumscribe(complex_exponential[6]))
        self.wait()

        r = 2.5 #circle radius
        ghost_circle = Circle(radius=r).set_opacity(0)
        circular_crystal = VGroup()
        N = 12 #Number of points on circle
        for i in range(N):
            angle = 2*PI/N * i
            circular_crystal.add(
                Dot(color=BLUE).move_to(ghost_circle.get_center() + r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        temp = self.mobjects
        self.play(
            ReplacementTransform(lattice_points_1D.copy(), circular_crystal),
            *[FadeOut(mob) for mob in self.mobjects]
        )
        self.wait()
        self.play(*[FadeIn(mob) for mob in temp], FadeOut(circular_crystal))
        self.wait()

        self.play(
            wave_func.animate.shift(DOWN*2.2),
            lattice_points_1D.animate.shift(DOWN*2.2),
            Unwrite(lattice_constant)
        )
        self.wait()
        self.play(
            ReplacementTransform(complex_exponential, nth_roots[1:3].shift(UP*0.5))
        )
        self.wait()
        self.play(
            Write(nth_roots[3].shift(UP*0.5))
        )
        self.wait(2)
        bloch_theorem = MathTex("\\psi(x+sa) = e^{","i2\\pi s/N","}\\psi(x)").move_to([0,-0.5,0])
        self.play(Write(bloch_theorem))
        self.wait(2)
        bloch_title = Tex("Bloch's Theorem").next_to(bloch_theorem, DOWN, buff=0.3)
        bloch_theorem = VGroup(bloch_theorem, bloch_title)
        bloch_rectangle = SurroundingRectangle(bloch_theorem, color=GOLD)
        bloch_theorem.add(bloch_rectangle)
        self.play(FadeIn(bloch_title), FadeIn(bloch_rectangle))
        self.wait()

        #Adding 'k'
        k = MathTex("^{ika}")
        k_expression = MathTex("k =", "2\\pi s/Na")
        bloch_theorem.add(k, k_expression)
        self.play(
            ReplacementTransform(bloch_theorem[0][1], k.move_to(bloch_theorem[0][1].get_center()+LEFT*0.3)),
            bloch_theorem[0][2].animate.shift(LEFT*0.6)
        )
        self.wait()
        new_rectangle = SurroundingRectangle(bloch_theorem, color=GOLD).move_to(bloch_title.get_center()+LEFT*0.2+UP*0.2)
        self.play(
            Write(k_expression.move_to(bloch_title.get_center()+UP*0.2)), 
            bloch_title.animate.shift(DOWN*0.3),
            ReplacementTransform(bloch_rectangle, new_rectangle)
        )
        self.wait()
        bloch_1D = Tex("Bloch's Theorem", " (In 1D)").move_to(bloch_title.get_center()+LEFT*0.2)
        self.play(
            bloch_title.animate.move_to(bloch_1D[0].get_center()),
            FadeIn(bloch_1D[1])
        )
        bloch_theorem.add(bloch_1D)

        self.wait(3)
        no_fade_group = VGroup(
            bloch_theorem, wave_func, lattice_points_1D, lattice_bracket
        )
        fade_group = Group()
        for mob in self.mobjects:
            if mob in no_fade_group or mob in bloch_theorem or mob is new_lb or mob in bloch_1D or mob is new_rectangle:
                continue
            fade_group.add(mob)
        self.play(
            *[
                Unwrite(mob) for mob in VGroup(nth_roots, complex_exponential, PSI[1],
                                               PSI[3], PSI[5], PSI[7], BVK)
            ],
            run_time = 0.7
        )
        self.wait()
        #Final form of bloch's theorem in 1D
        satisfactoy_wf = MathTex(
            "\\psi(x)",     #0
            "=",            #1
            "e^{ikx}",      #2
            "\\implies",    #3
            "\\psi(x+",     #4
            "s",            #5
            "a)",           #6
            "=",            #7
            "e^{i2\\pi ",   #8
            "s^2",            #9
            "/N}",          #10             
            "e^{ikx}",      #11   
        ).shift(UP*2)
        psi_replacement = MathTex("N", "^{sN}", "\\psi(x)")
        self.play(Write(satisfactoy_wf[0:3]))
        self.wait()
        self.play(Write(satisfactoy_wf[3:]))
        self.wait()
        self.play(
            ReplacementTransform(satisfactoy_wf[5], psi_replacement[0].move_to(satisfactoy_wf[5].get_center()+UP*0.05+LEFT*0.03)),
            ReplacementTransform(satisfactoy_wf[9], psi_replacement[1].move_to(satisfactoy_wf[9].get_center()+DOWN*0.05+RIGHT*0.1)),
            satisfactoy_wf[10:].animate.shift(RIGHT*0.1)
        )
        self.wait()
        self.play(
            ReplacementTransform(satisfactoy_wf[11], psi_replacement[2].move_to(satisfactoy_wf[11].get_center()+RIGHT*0.25)),
        )
        self.wait()
        self.play(
            FadeOut(satisfactoy_wf[8]), FadeOut(satisfactoy_wf[10]), FadeOut(psi_replacement[1]),
            psi_replacement[2].animate.shift(LEFT*1.55+DOWN*0.05)
        )
        self.wait()

class ElectronWaveFunction(Scene):
    def construct(self):
        lattice_points_1D = VGroup()
        ah_lattice = 1 #Lattice constant horizontal
        point_count_h = 20 #Horizontal point count

        for i in range(point_count_h):
            lattice_points_1D.add(Dot(color=BLUE).move_to(LEFT*(10+2*ah_lattice)+RIGHT*ah_lattice*i))
        lattice_points_1D = lattice_points_1D.move_to([0,0,0])

        self.play(FadeIn(lattice_points_1D))
        self.wait()

        lattice_bracket = always_redraw(
            lambda: BraceBetweenPoints(
                lattice_points_1D[int(len(lattice_points_1D)/2)-1].get_center(),
                lattice_points_1D[int(len(lattice_points_1D)/2)].get_center()
                ).set_color(RED)
        )

        self.play(Write(lattice_bracket))
        self.wait()

        lattice_constant_val = ValueTracker(np.linalg.norm(lattice_points_1D[0].get_center() - lattice_points_1D[1].get_center())) 

        lattice_constant = VGroup(
            MathTex("a"),
            DecimalNumber(lattice_constant_val.get_value(), num_decimal_places=3, include_sign=False, unit=None)
        ).move_to([0,.8,0]).set_color(RED).scale(1.5)
        lattice_constant[1] = lattice_constant[1].next_to(lattice_constant[0], DOWN, buff=0.2)        
        lattice_constant[1].add_updater(
            lambda d: d.set_value(lattice_constant_val.get_value())
        )
        lattice_constant.add(SurroundingRectangle(lattice_constant, color=DARK_BLUE))
        
        self.play(lattice_points_1D.animate.shift(DOWN*1))
        self.wait()
        self.play(Write(lattice_constant))
        self.wait()
        self.play(lattice_points_1D.animate.scale(1.5), lattice_constant_val.animate.set_value(1.5))
        self.wait()
        self.play(lattice_points_1D.animate.scale(0.8), lattice_constant_val.animate.set_value(1.5*0.8))
        self.wait()
        new_lb = always_redraw(lambda: lattice_bracket.copy().shift(RIGHT*lattice_constant_val.get_value()))
        self.play(Write(new_lb))
        self.wait()
        self.play(lattice_points_1D.animate.shift(LEFT*lattice_constant_val.get_value()))
        self.remove(lattice_bracket)
        self.wait()
        self.play(
            lattice_constant.animate.next_to(new_lb.get_tip(), DOWN),
        )

        #Adding sine wave
        wavelength = 2
        ax = Axes(
            x_range = [-len(lattice_points_1D)/2+1, len(lattice_points_1D)/2],
            y_range = [-2,2],
            x_length = lattice_points_1D[-1].get_center()[0] - lattice_points_1D[0].get_center()[0],
            y_length = 4
        ).shift(DOWN*1)
        wave_func = ax.plot(lambda x: np.cos(2*PI/wavelength * x)**2, color=BLUE_E)

        self.play(Create(wave_func), run_time=2)

        self.wait()
        self.play(wave_func.animate.shift(RIGHT*lattice_constant_val.get_value()), run_time=2)
        self.wait()

        BVK = MathTex("\\textbf{Translational Symmetry}").shift(UP*1.5).scale(1.2)

        self.play(Write(BVK))
        self.wait()
        self.play(
            BVK.animate.scale(0.8).move_to([0,3,0]),
        )

        PSI = MathTex(
            "|",                                #0
            "\\psi(x)",                            #1
            "|^2",                              #2
            "=",                       #3     
            "|",                         #4
            "\\psi(x+sa)",                #5
            "|^2",                     #6
            ", s = 0, 1, 2... ",        #7
            # "\\text{ Number of atoms in crystal}" #8
        ).move_to([0,1.5,0]).scale(0.9)

        statement = Text("Electron Wave Function").scale(0.8).next_to(PSI[1], UP, buff=0.2)
        
        self.play(Write(PSI[1]), Write(statement))
        self.wait()
        self.play(FadeOut(statement))
        self.wait()
