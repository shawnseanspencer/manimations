from manim import *
from manim_physics import *

class MultiPendulumExample(SpaceScene):
    def construct(self):
        p = MultiPendulum(RIGHT, LEFT+UP*1.5)
        self.add(p)
        self.make_rigid_body(*p.bobs)
        p.start_swinging()
        self.add(TracedPath(p.bobs[-1].get_center, stroke_color=BLUE))
        self.wait(10)

def get_oscillator(
    start = [0,0,0],
    length: float = 4,
    axis: np.array = np.array([1, 0, 0]),
    k: float = 5,
    block_side_length: float = 1,
    color = WHITE,
    spring_width: float = 0.7
):

    #Normalizing axis
    if np.linalg.norm(axis) != 1:
        axis = axis/np.linalg.norm(axis)

    angle = np.arctan(axis[1]/axis[0])

    temp_ax = Axes(
        x_range = [0, 2*PI],
        y_range = [-1, 1],
        x_length = length,
        y_length = spring_width
    ).rotate(angle).set_opacity(0)
     
    block = Square(
        side_length = block_side_length
    ).set_fill(color).move_to([0,0,0])

    temp_line = always_redraw(
        lambda: Line(
            start = start,
            end = block.get_left(),
        ).set_opacity(0)
    )

    temp_ax = always_redraw(
            lambda: NumberPlane(
            x_range = [0, 2*PI],
            y_range = [-1, 1],
            x_length = np.linalg.norm(temp_line.get_left() - temp_line.get_right()),
            y_length = spring_width
        ).move_to(temp_line.get_center()).set_opacity(0)
    )

    spring_period =  temp_ax.p2c(block.get_left())[0] / 5
    spring = temp_ax.plot(
        lambda x: np.sin(2*PI/spring_period * x), 
        x_range = [0, temp_ax.p2c(block.get_left())[0]],
    ).set_color(color)
    spring.add_updater(
            lambda mob: mob.become(temp_ax.plot(
            lambda x: np.sin(2*PI/spring_period * x), 
            x_range = [0, temp_ax.p2c(block.get_left())[0]],
        ).set_color(color))
    )

    return {
        "system": VGroup(block, spring, temp_line, temp_ax),
        "attrs": [length, k]
    }


class EmphasisWords(Scene):
    def construct(self):
        language = VGroup(
            Text("What is the", font = "PT Mono"),
            Text("language", font = "PT Mono"),
            Text("of Nature?", font = "PT Mono")
        ).arrange(DOWN).scale(2)
        for item in language:
            self.play(AddTextLetterByLetter(item))
        self.wait()
        self.play(FadeOut(language))


        DE = Text("Differential \n Equations", font = "Chalkduster").scale(2)
        self.play(AddTextLetterByLetter(DE), run_time = 2)
        self.wait()
        self.play(DE.animate.shift(UP*1.5))
        DE_ex = MathTex(
            "F = ma = m\\frac{d^2x}{dt^2} = \\frac{dp}{dt}"
        ).scale(2).next_to(DE,DOWN)
        self.play(Write(DE_ex))
        self.wait()
        self.play(FadeOut(DE), FadeOut(DE_ex))

        unsolvable = Text("Unsolvable", font = "Chalkduster").scale(2)
        self.play(AddTextLetterByLetter(unsolvable), run_time = 2)
        self.wait()
        self.play(unsolvable.animate.shift(UP*1.5))
        exception = Text("but... we can approximate").next_to(unsolvable, DOWN, buff = 0.8)
        self.play(Write(exception))
        self.wait()
        self.play(FadeOut(unsolvable), FadeOut(exception))
        self.wait()


class Oscillator(Scene):
    def construct(self):
        setup = get_oscillator(start = [-3, 0, 0], length = 5)
        system = setup['system']

        self.add(system)
        self.wait()

        time = np.linspace(-4*PI, 4*PI, 300)

        dt = 0.02
        for t in time:
            self.play(system[0].animate.set_x(2*np.sin(t)), rate_func = linear, run_time = dt)

        self.wait()

class Test(Scene):
    def construct(self):
        time = np.linspace(-10, 10, 100)
        cube = Square()
        dt = 0.1

        def pos_func(t):
            return np.cos(t)

        for t in time:
            self.play(cube.animate.move_to([pos_func(t), 0, 0]), rate_func = linear, run_time = dt)
        self.wait()

class Spring(Scene):
    def construct(self):

        ax = Axes(
            x_range = [-PI, 4],
            x_length = 5,
            y_length = 3
        ).shift(LEFT*2)
        func = always_redraw(lambda: ax.plot(
            lambda x: np.sin(6*x),
            x_range = [-PI, PI]
        ))
        block = always_redraw(
            lambda: Cube(side_length = 1).move_to(ax.c2p(PI,0,0)+0.5*RIGHT)
        )
        self.add(block)

        self.add(func)
        self.play(ax.animate.stretch_about_point(2, 0, [ax.c2p(-PI,0,0)]))
        self.wait()

class BlockSpring(VGroup):
    def __init__(self,
                start: np.array = np.array([0,0,0]),
                length: float = 4,
                axis: np.array = np.array([1, 0, 0]),
                k: float = 5,
                block_side_length: float = 1,
                color = WHITE,
                spring_width: float = 0.7,
                spring_period: float = 2,
                **kwargs
        ):
        self.length = length
        self.startP = start
        self.spring_width = spring_width
        self.spring_period = spring_period
        self.axis = axis/np.linalg.norm(axis)
        self.CustColor = color
        self.block = Square(side_length=block_side_length, color=self.CustColor)
        self.block.shift((length-0.5*block_side_length)*RIGHT)
        self.block.rotate(angle=np.arctan2(axis[1],axis[0]),about_point=ORIGIN)
        self.block.shift(start)
        self.spring = VMobject()
        super().__init__(self.block, self.spring, **kwargs)
        self.add_updater(self.springUpdater, call_updater=True)

    def springUpdater(self, mobj):
        verts = self.block.get_vertices()
        endP = 0.5*(verts[1]+verts[2])
        helpline = Line(self.startP, endP)
        ax = Axes(
            x_range = [0, 2*PI],
            y_range = [-1, 1],
            x_length = np.linalg.norm(endP-self.startP),
            y_length = self.spring_width
        )
        ax.shift(self.startP-ax.c2p(0,0))
        ax.rotate(helpline.get_angle(), about_point=ax.c2p(0,0))
        self.spring.become(ax.plot(lambda x: np.sin(self.spring_period*x),color=self.CustColor))
    def move_block_to(self, point):
        verts = self.block.get_vertices()
        endP = 0.5*(verts[1]+verts[2])
        self.block.shift(point-endP)
    def set_length(self, length):
        verts = self.block.get_vertices()
        endP = 0.5*(verts[1]+verts[2])
        point = self.startP + length*self.axis
        self.block.shift(point-endP)

class Thumbnail(Scene):
    def construct(self):
        bs1 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=WHITE, spring_width = 1.3, block_side_length = 2, spring_period=3, length = 8)
        self.add(bs1)

class testBlockSpring(Scene):
    def construct(self):
        bs1 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs1[0].shift(RIGHT)

        def pos_func(t):
            return 2*np.cos(t)

        self.play(Create(bs1))
        self.wait()
        time = 4*PI
        dt = 0.1
        t = 0
        while t < time:
            self.play(
                bs1[0].animate.move_to([pos_func(t), 0, 0]), 
                rate_func = linear, 
                run_time = dt/2
            )
            t += dt
        self.wait()

        harmonic_oscillator = Tex("Harmonic Oscillator").scale(1.5).shift(UP*2)
        self.play(Write(harmonic_oscillator))
        self.wait()

        bs2 = BlockSpring(start = [-3.5,-2.5,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs2[0].shift(RIGHT)
        self.play(ReplacementTransform(bs1, bs2))
        self.wait()

        #Hooke's Law
        hooke = MathTex("-kx = ma \\implies \\frac{d^2x}{dt^2} = -\\frac{k}{m} x").scale(1.5)
        self.play(Write(hooke))
        self.wait()
        self.play(Unwrite(hooke), Unwrite(harmonic_oscillator))
        self.wait()

        upward_mover = ValueTracker(bs2[0].get_center()[1])
        path_head = always_redraw(
            lambda: Dot().move_to([bs2[0].get_center()[0], upward_mover.get_value(), 9])
        )
        path_tracer = TracedPath(path_head.get_center, dissipating_time=100, stroke_opacity = 0.38, stroke_width=8).set_color(GREEN)
        connecting_line = always_redraw(
            lambda: DashedLine(
                start = path_head.get_center(),
                end = bs2[0].get_center()
            ).set_opacity(0.5)
        )
        self.add(connecting_line)
        self.add(path_head)
        self.add(path_tracer)

        t = 0
        while t < time + 2*PI:
            self.play(
                bs2[0].animate.move_to([pos_func(t), -2.5, 0]), 
                upward_mover.animate.set_value(bs2[0].get_center()[1] + t/2.5),
                rate_func = linear, 
                run_time = dt/2
            )
            t += dt
        self.wait()

        #Sinusoidal Solution
        sinusoid = VGroup(
            MathTex(
                "x(t)",                 #0
                "=",                    #1
                "A",                    #2
                "\\cos",                #3  
                "\\left(",              #4
                "\\sqrt{",              #5
                "k",                    #6
                "\\over",               #7
                "m}",                   #8
                "t",                    #9
                "-",                    #10
                "\\phi",                #11
                "\\right)"              #12
            ),
            MathTex("v(t) = -A\\sqrt{\\frac{k}{m}}\\cos\\left(\\sqrt{\\frac{k}{m}}t - \\phi\\right)")
        ).arrange(DOWN).shift(UP*0.5)

        sinusoid = VGroup(BackgroundRectangle(sinusoid), *[mob for mob in sinusoid])
        self.play(Write(sinusoid))
        self.wait()

        self.play(Create(A_rect := SurroundingRectangle(sinusoid[1][2], color = PURPLE)))
        self.wait()
        self.play(Create(m_rect := SurroundingRectangle(sinusoid[1][8], color = RED)))
        self.wait()
        self.play(Create(k_rect := SurroundingRectangle(sinusoid[1][6], color = BLUE)))
        self.wait()
        self.play(FadeOut(VGroup(A_rect, m_rect, k_rect)))
        self.wait()
        strikeout_lines = VGroup(
            Line(start=sinusoid[1].get_left(), end=sinusoid[1].get_right(), stroke_width = 3),
            Line(start=sinusoid[2].get_left(), end=sinusoid[2].get_right(), stroke_width = 3)
        )
        self.play(Create(strikeout_lines))
        self.wait()

        self.play(
            Unwrite(
                VGroup(sinusoid, strikeout_lines, connecting_line, path_tracer, path_head)
            )
        )
        self.wait()


        """
        self.play(
            bs1[0].animate.shift(3*bs1.axis),
            bs2.animate.set_length(6),
            bs3.animate.move_block_to(np.array([-6,-1,0])),
            run_time=3
        )
        self.play(
            bs1[0].animate.shift(-5*bs1.axis),
            bs2.animate.set_length(1),
            bs3.animate.move_block_to(np.array([0,-1,0])),
            run_time=3
        )
        """

class EulerDerivation(Scene):
    def construct(self):
        derivation = BulletedList(
            "$$F = ma \\implies -kx = ma \\implies a = -\\frac{kx}{m}$$",
            "Split up time into discrete steps $\\Delta t$",
            "Iteratively calculate the next $v$ and $x$",
            "$$v_{i+1} = v_i + a\\Delta t$$",
            "$$x_{i+1} = x_i + v_{i}\\Delta t$$"
        )

        for item in derivation:
            self.play(Write(item))
            self.wait()

        self.play(*[FadeOut(mob) for mob in self.mobjects])
        self.wait()

        generic = VGroup(
            MathTex("\\frac{dx}{dt} = f(x, t), \\text{ } x(t_0) = x_0"),
            Tex("Euler's Method"),
            MathTex("x_{i+1} = x_i + f(x_i, t_i)\\Delta t"),
            Tex("Implicit Euler's Method"),
            MathTex("x_{i+1}", "=", "x_i","+", "f(x_{i+1}, t_{i+1})\\Delta t"),
        ).arrange(DOWN, buff = 0.5)
        for i, item in enumerate(generic):
            self.play(Write(item))
            self.wait()
            if i == 1:
                item.add(Tex("(Initial Value Problems)").next_to(item, RIGHT))
                self.play(Write(item[-1]))
                self.play(item.animate.set_x(0))
                self.wait()
            elif i == 4:

                self.play(item[0].animate.set_x(item[2].copy().get_center()[0]), item[2].animate.set_x(item[0].copy().get_center()[0]), item[4].animate.shift(
                    RIGHT * np.linalg.norm((item[0].get_right()-item[0].get_left()) - (item[2].get_right() - item[2].get_left()))/2
                ), item[1].animate.shift(LEFT*np.linalg.norm((item[0].get_right()-item[0].get_left()) - (item[2].get_right() - item[2].get_left()))/2),
                ReplacementTransform(item[3], temp_minus := MathTex("-").move_to(item[3].copy().get_center()+RIGHT * np.linalg.norm((item[0].get_right()-item[0].get_left()) - (item[2].get_right() - item[2].get_left()))/2))
                )
                generic.add(temp_minus)
                self.wait()

                generic[3].add(Tex("(Backwards)").next_to(generic[3], RIGHT))
                self.play(Write(generic[3][-1]))
                self.play(generic[3].animate.set_x(0))
                self.wait()

class RKDerivation(Scene):
    def construct(self):
        method = BulletedList(
            "Average of the two methods,... well almost",
            "The energy increase and loss might cancel out?",
            "$$k_1 = f(x_i, t_i), \\text{ } k_2 = f(x_i + k_1, t_i + \\Delta t)$$",
            "$x_{i} + k_1 = x_i + f(x_i, t_i)\\Delta t$ same as in Euler's Method",
            "$$\\text{But now, let }x_{i+1} = x_i +\\frac{dt}{2 }(k_1 + k_2)$$"
        ).scale(0.9)
        for i, item in enumerate(method):
            self.play(Write(item))
            self.wait()


class EulerDemonstration(Scene):
    def construct(self):

        Amp = 2
        dt = 0.1
        k = 5
        m = 1

        ax = Axes(
            x_range = [0, 2*PI],
            y_range = [-2, 2],
            x_length = 7,
            y_length = 3
        )
        labels = ax.get_axis_labels(x_label="t", y_label = "x")

        def spring_func(t):
            return Amp*np.cos(np.sqrt(k/m) * t)
        def spring_deriv(t):
            return -Amp*np.sqrt(k/m)*np.sin(np.sqrt(k/m) * t)
        def accel(x):
            return -k * x / m
    
        true_func = ax.plot(
            lambda t: spring_func(t),
            color = RED
        )

        self.play(Create(VGroup(ax, labels)), run_time = 2)
        self.play(Create(true_func), run_time = 1.5)

        tracker_point = Dot(color = BLUE).move_to(ax.c2p(0,2,0)).scale(0.9)
        v_line = Line(
            start = tracker_point.get_center(),
            end = ax.c2p(
                0 + dt, spring_func(0) + (spring_deriv(0)) * dt
            )
        ).set_opacity(0.5)

        self.play(
            VGroup(ax, labels, true_func).animate.scale(5, about_point = tracker_point.get_center()),
            run_time = 2
        )
        self.wait()

        self.play(Create(VGroup(tracker_point.scale(5), v_line.scale(5, about_point = tracker_point.get_center()))))
        self.wait()
        self.play(tracker_point.animate.move_to(v_line.get_end()))
        self.wait()


        self.play(
            VGroup(ax, labels, true_func, tracker_point, v_line).animate.scale(1/5, about_point = tracker_point.get_center()),
            run_time = 2
        )
        self.wait()

        approx_lines = VGroup(v_line)
        tracker_points = VGroup(tracker_point, tracker_point.copy())

        def euler_step(x, v, dt, k, m):
            a = - (k / m) * x
            v_new = v + a * dt
            x_new = x + v * dt
            return x_new, v_new

        final_time = 6
        t = dt
        x_euler, v_euler = Amp, 0
        x_euler, v_euler = euler_step(x_euler, v_euler, dt, k, m)
        while t < final_time:
            
            x_euler, v_euler = euler_step(x_euler, v_euler, dt, k, m)

            approx_lines.add(
                Line(
                    start = tracker_point.get_center(),
                    end = ax.c2p(
                        t + dt, 
                        x_euler
                    )
                ).set_opacity(0.5)
            )
            self.play(Create(approx_lines[-1]), run_time = 0.1)
            self.play(tracker_point.animate.move_to(approx_lines[-1].get_end()), run_time = 0.1)
            tracker_points.add(tracker_point.copy())
            self.add(tracker_points[-1])
            t += dt

        self.wait()
        self.play(VGroup(tracker_points, approx_lines).animate.set_opacity(0.3))
        self.wait()
        
        better_lines = VGroup()
        better_point = Dot(color = PURPLE).move_to(ax.c2p(0,2,0)).scale(0.7)
        better_points = VGroup()
        small_dt = 0.02
        t = 0
        x_euler, v_euler = Amp, 0
        while t < final_time:
            x_euler, v_euler = euler_step(x_euler, v_euler, small_dt, k, m)

            better_lines.add(
                Line(
                    start = better_point.get_center(),
                    end = ax.c2p(
                        t + small_dt, 
                        x_euler
                    )
                ).set_opacity(0.5).set_color(PURPLE)
            )
            self.play(Create(better_lines[-1]), run_time = 0.03)
            self.play(better_point.animate.move_to(better_lines[-1].get_end()), run_time = 0.1)
            better_points.add(better_point.copy())
            self.add(better_points[-1])
            t += small_dt
        self.wait()


        worse_lines = VGroup()
        worse_point = Dot(color = GREEN).move_to(ax.c2p(0,2,0)).scale(0.9)
        worse_points = VGroup()
        big_dt = 0.3
        t = 0

        x_euler, v_euler = Amp, 0
        while t < final_time:
            x_euler, v_euler = euler_step(x_euler, v_euler, big_dt, k, m)

            worse_lines.add(
                Line(
                    start = worse_point.get_center(),
                    end = ax.c2p(
                        t + big_dt, 
                        x_euler
                    )
                ).set_opacity(0.5).set_color(GREEN)
            )
            self.play(Create(worse_lines[-1]), run_time = 0.1)
            self.play(worse_point.animate.move_to(worse_lines[-1].get_end()), run_time = 0.1)
            worse_points.add(worse_point.copy())
            self.add(worse_points[-1])
            t += big_dt

        self.wait()

        time_steps = MathTex("\\Delta t = 0.1,", "\\text{ }","\\Delta t = 0.3", "\\text{ }", "\\Delta t = 0.02").shift(UP*2.5).scale(1.2)
        time_steps[0].set_color(BLUE)
        time_steps[2].set_color(GREEN)
        time_steps[4].set_color(PURPLE)
        self.play(Write(time_steps))
        self.wait()


class ImplicitDemonstration(Scene):
    def construct(self):
        Amp = 2
        dt = 0.1
        k = 5
        m = 1

        ax = Axes(
            x_range = [0, 2*PI],
            y_range = [-2, 2],
            x_length = 7,
            y_length = 3
        )
        labels = ax.get_axis_labels(x_label="t", y_label = "x")

        def spring_func(t):
            return Amp*np.cos(np.sqrt(k/m) * t)
    
        true_func = ax.plot(
            lambda t: spring_func(t),
            color = RED
        )

        self.play(Create(VGroup(ax, labels)), run_time = 2)
        self.play(Create(true_func), run_time = 1.5)

        tracker_point = Dot(color = YELLOW).move_to(ax.c2p(0,2,0)).scale(0.9)
        self.wait()

        approx_lines = VGroup()
        tracker_points = VGroup(tracker_point, tracker_point.copy())

        def euler_step(x, v, dt, k, m):
            a = - (k / m) * x
            v_new = v + a * dt
            x_new = x + v * dt
            return x_new, v_new
        
        def implicit_step(x, v, dt, k, m):
            A = 1 + (dt**2 * k / m)
            B = -dt * k / m

            x_new = (x + dt * v) / A
            v_new = v + dt * (-k / m * x_new)
            
            return x_new, v_new

        final_time = 6
        t = dt
        x_euler, v_euler = Amp, 0
        x_euler, v_euler = implicit_step(x_euler, v_euler, dt, k, m)
        while t < final_time:
            
            x_euler, v_euler = implicit_step(x_euler, v_euler, dt, k, m)

            approx_lines.add(
                Line(
                    start = tracker_point.get_center(),
                    end = ax.c2p(
                        t + dt, 
                        x_euler
                    )
                ).set_opacity(0.5)
            )
            self.play(Create(approx_lines[-1]), run_time = 0.1)
            self.play(tracker_point.animate.move_to(approx_lines[-1].get_end()), run_time = 0.1)
            tracker_points.add(tracker_point.copy())
            self.add(tracker_points[-1])
            t += dt

class RKDemonstration(Scene):
    def construct(self):
        Amp = 2
        dt = 0.1
        k = 5
        m = 1

        ax = Axes(
            x_range = [0, 2*PI],
            y_range = [-2, 2],
            x_length = 7,
            y_length = 3
        )
        labels = ax.get_axis_labels(x_label="t", y_label = "x")

        def spring_func(t):
            return Amp*np.cos(np.sqrt(k/m) * t)
    
        true_func = ax.plot(
            lambda t: spring_func(t),
            color = RED
        )

        self.play(Create(VGroup(ax, labels)), run_time = 2)
        self.play(Create(true_func), run_time = 1.5)

        tracker_point = Dot(color = ORANGE).move_to(ax.c2p(0,2,0)).scale(0.9)
        self.wait()

        approx_lines = VGroup()
        tracker_points = VGroup(tracker_point, tracker_point.copy())

        def euler_step(x, v, dt, k, m):
            a = - (k / m) * x
            v_new = v + a * dt
            x_new = x + v * dt
            return x_new, v_new
        
        def implicit_step(x, v, dt, k, m):
            A = 1 + (dt**2 * k / m)
            B = -dt * k / m

            x_new = (x + dt * v) / A
            v_new = v + dt * (-k / m * x_new)
            
            return x_new, v_new
        
        def rk2_step(x, v, dt, k, m):
            def derivatives(x, v):
                return v, - (k / m) * x
            
            k1x, k1v = derivatives(x, v)
            k2x, k2v = derivatives(x + dt*k1x, v + dt*k1v)
            
            x_new = x + (dt / 2) * (k1x + k2x)
            v_new = v + (dt / 2) * (k1v + k2v)
    
            return x_new, v_new
        
        def rk4_step(x, v, dt, k, m):
            def derivatives(x, v):
                return v, - (k / m) * x
            
            k1x, k1v = derivatives(x, v)
            k2x, k2v = derivatives(x + 0.5 * dt * k1x, v + 0.5 * dt * k1v)
            k3x, k3v = derivatives(x + 0.5 * dt * k2x, v + 0.5 * dt * k2v)
            k4x, k4v = derivatives(x + dt * k3x, v + dt * k3v)
            
            x_new = x + (dt / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
            v_new = v + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
            
            return x_new, v_new

        final_time = 6
        t = dt
        x_euler, v_euler = Amp, 0
        x_euler, v_euler = rk2_step(x_euler, v_euler, dt, k, m)
        while t < final_time:
            
            x_euler, v_euler = rk2_step(x_euler, v_euler, dt, k, m)

            approx_lines.add(
                Line(
                    start = tracker_point.get_center(),
                    end = ax.c2p(
                        t + dt, 
                        x_euler
                    )
                ).set_opacity(0.5)
            )
            self.play(Create(approx_lines[-1]), run_time = 0.1)
            self.play(tracker_point.animate.move_to(approx_lines[-1].get_end()), run_time = 0.1)
            tracker_points.add(tracker_point.copy())
            self.add(tracker_points[-1])
            t += dt
        self.wait()

        rk4_point = Dot(color = GOLD).move_to(ax.c2p(0,2,0)).scale(0.9)
        self.wait()

        self.play(VGroup(tracker_points, approx_lines).animate.set_opacity(0.3))
        self.wait()

        rk4_lines = VGroup()
        rk4_points = VGroup(rk4_point, rk4_point.copy())

        t = dt
        x_euler, v_euler = Amp, 0
        x_euler, v_euler = rk4_step(x_euler, v_euler, dt, k, m)
        while t < final_time:
            
            x_euler, v_euler = rk4_step(x_euler, v_euler, dt, k, m)

            rk4_lines.add(
                Line(
                    start = rk4_point.get_center(),
                    end = ax.c2p(
                        t + dt, 
                        x_euler
                    )
                ).set_opacity(0.5)
            )
            self.play(Create(rk4_lines[-1]), run_time = 0.1)
            self.play(rk4_point.animate.move_to(rk4_lines[-1].get_end()), run_time = 0.1)
            rk4_points.add(rk4_point.copy())
            self.add(rk4_points[-1])
            t += dt

class EulerSpring(Scene):
    def construct(self):

        Amp = 2
        dt = 0.05
        k = 5
        m = 1

        bs1 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs1[0].shift(RIGHT)

        def pos_func(t):
            return 2*np.cos(t)

        self.play(Create(bs1))
        self.wait()
        time = 4*PI

        bs2 = BlockSpring(start = [-3.5,-2.5,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs2[0].shift(RIGHT)
        self.play(ReplacementTransform(bs1, bs2))
        self.wait()

        #Hooke's Law
        upward_mover = ValueTracker(bs2[0].get_center()[1])
        path_head = always_redraw(
            lambda: Dot().move_to([bs2[0].get_center()[0], upward_mover.get_value(), 9])
        )
        path_tracer = TracedPath(path_head.get_center, dissipating_time=100, stroke_opacity = 0.38, stroke_width=8).set_color(GREEN)
        connecting_line = always_redraw(
            lambda: DashedLine(
                start = path_head.get_center(),
                end = bs2[0].get_center()
            ).set_opacity(0.5)
        )
        self.add(connecting_line)
        self.add(path_head)
        self.add(path_tracer)

        def euler_step(x, v, dt, k, m):
            a = - (k / m) * x
            v_new = v + a * dt
            x_new = x + v * dt
            return x_new, v_new

        x_euler, v_euler = Amp, 0
        t = 0
        while t < time + 2*PI:

            x_euler, v_euler = euler_step(x_euler, v_euler, dt, k, m)

            self.play(
                bs2[0].animate.move_to([x_euler, -2.5, 0]), 
                upward_mover.animate.set_value(bs2[0].get_center()[1] + t/2.5),
                rate_func = linear, 
                run_time = dt/2
            )
            t += dt
        self.wait()

class ImplicitSpring(Scene):
    def construct(self):

        Amp = 2
        dt = 0.05
        k = 5
        m = 1

        bs1 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs1[0].shift(RIGHT)

        def pos_func(t):
            return 2*np.cos(t)

        self.play(Create(bs1))
        self.wait()
        time = 4*PI

        bs2 = BlockSpring(start = [-3.5,-2.5,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs2[0].shift(RIGHT)
        self.play(ReplacementTransform(bs1, bs2))
        self.wait()

        #Hooke's Law
        upward_mover = ValueTracker(bs2[0].get_center()[1])
        path_head = always_redraw(
            lambda: Dot().move_to([bs2[0].get_center()[0], upward_mover.get_value(), 9])
        )
        path_tracer = TracedPath(path_head.get_center, dissipating_time=100, stroke_opacity = 0.38, stroke_width=8).set_color(GREEN)
        connecting_line = always_redraw(
            lambda: DashedLine(
                start = path_head.get_center(),
                end = bs2[0].get_center()
            ).set_opacity(0.5)
        )
        self.add(connecting_line)
        self.add(path_head)
        self.add(path_tracer)

        def euler_step(x, v, dt, k, m):
            a = - (k / m) * x
            v_new = v + a * dt
            x_new = x + v * dt
            return x_new, v_new
        
        def implicit_step(x, v, dt, k, m):
            A = 1 + (dt**2 * k / m)
            B = -dt * k / m

            x_new = (x + dt * v) / A
            v_new = v + dt * (-k / m * x_new)
            
            return x_new, v_new

        x_euler, v_euler = Amp, 0
        t = 0
        while t < time + 2*PI:

            x_euler, v_euler = implicit_step(x_euler, v_euler, dt, k, m)

            self.play(
                bs2[0].animate.move_to([x_euler, -2.5, 0]), 
                upward_mover.animate.set_value(bs2[0].get_center()[1] + t/2.5),
                rate_func = linear, 
                run_time = dt/2
            )
            t += dt
        self.wait()

class ForceIncrementSpring(Scene):
    def construct(self):

        bs1 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs1[0].shift(RIGHT)
        ghost = bs1[0].copy().set_stroke(opacity=0.5)
        self.add(bs1)
        self.add(ghost)

        parameters = VGroup(
            MathTex("\\Delta t, k, m"),   #0,
            MathTex("F_0"),         #1,
            MathTex("v_0"),         #2,
            MathTex("x_0")          #3 
        ).arrange(RIGHT, buff = 0.3).shift(UP*2.5)
        parameters[1].set_color(RED)
        parameters[2].set_color(BLUE)
        parameters[3].set_color(PURPLE)

        self.play(Write(parameters))
        self.wait()

        F_arrow = Arrow(
            start = bs1[0].get_center() + UP*0.8 + RIGHT*0.2,
            end = bs1[0].get_center() + UP*0.8 + LEFT*1.2,
            stroke_width = 6,
            color = RED
        )
        F_arrow_label = MathTex("\\vec{F}").next_to(F_arrow, UP).set_color(RED)
        
        self.play(GrowArrow(F_arrow))
        self.play(Write(F_arrow_label))
        self.wait()

        v_arrow = Arrow(
            start = bs1[0].get_center() + DOWN*0.8 + RIGHT*0.2,
            end = bs1[0].get_center() + DOWN*0.8 + LEFT*1.1,
            stroke_width = 6,
            color = BLUE
        )
        v_arrow_label = MathTex("\\vec{v}").next_to(v_arrow, DOWN).set_color(BLUE)

        self.play(GrowArrow(v_arrow))
        self.play(Write(v_arrow_label))
        self.wait()

        v_value_label = MathTex("\\Delta\\vec v = \\frac{\\vec F}{m} \\Delta t").set_color(BLUE).move_to([4.5, 1, 0])
        v_value_label = VGroup(v_value_label, MathTex("\\left(\\vec{a} = \\frac{\\vec F}{m}\\right)").scale(0.7).next_to(v_value_label, UP))
        self.play(Write(v_value_label))
        self.wait()

        self.play(
            bs1[0].animate.shift(LEFT*1.1),
            rate_func = linear
        )
        self.wait()

        distance_line = DashedLine(
            start = bs1[0].get_center(),
            end = ghost.get_center(),
            color = PURPLE
        )
        distance_line_ghost = distance_line.copy().set_opacity(0.6)
        self.play(Create(distance_line))
        self.add(distance_line_ghost)
        self.wait()

        self.play(distance_line.animate.move_to([4.5, -0.3, 0]))
        self.wait()

        distance_brace = BraceLabel(
            distance_line,
            "\\Delta\\vec x = \\vec v \\Delta t",
            brace_direction = DOWN,
        ).set_color(PURPLE)
        self.play(Create(distance_brace))
        self.wait()

class ForceIncrementSpringEnding(Scene):
    def construct(self):

        bs1 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=WHITE, spring_width = 0.5, spring_period=3, length = 5)
        bs1[0].shift(RIGHT)
        ghost = bs1[0].copy().set_stroke(opacity=0.5)
        self.add(bs1)
        self.add(ghost)

        parameters = VGroup(
            MathTex("\\Delta t, k, m"),   #0,
            MathTex("F_0"),         #1,
            MathTex("v_0"),         #2,
            MathTex("x_0")          #3 
        ).arrange(RIGHT, buff = 0.3).shift(UP*2.5)
        parameters[1].set_color(RED)
        parameters[2].set_color(BLUE)
        parameters[3].set_color(PURPLE)

        self.play(Write(parameters))
        self.wait()

        F_arrow = Arrow(
            start = bs1[0].get_center() + UP*0.8 + RIGHT*0.2,
            end = bs1[0].get_center() + UP*0.8 + LEFT*1.2,
            stroke_width = 6,
            color = RED
        )
        F_arrow_label = MathTex("\\vec{F}").next_to(F_arrow, UP).set_color(RED)
        
        self.play(GrowArrow(F_arrow))
        self.play(Write(F_arrow_label))
        self.wait()

        v_arrow = Arrow(
            start = bs1[0].get_center() + DOWN*0.8 + RIGHT*0.2,
            end = bs1[0].get_center() + DOWN*0.8 + LEFT*1.1,
            stroke_width = 6,
            color = BLUE
        )
        v_arrow_label = MathTex("\\vec{v}").next_to(v_arrow, DOWN).set_color(BLUE)

        self.play(GrowArrow(v_arrow))
        self.play(Write(v_arrow_label))
        self.wait()

        v_value_label = MathTex("\\Delta\\vec v = \\frac{\\vec F}{m} \\Delta t").set_color(BLUE).move_to([4.5, 1, 0])
        v_value_label = VGroup(v_value_label, MathTex("\\left(\\vec{a} = \\frac{\\vec F}{m}\\right)").scale(0.7).next_to(v_value_label, UP))
        self.play(Write(v_value_label))
        self.wait()

        self.play(
            bs1[0].animate.shift(LEFT*1.1),
            rate_func = linear
        )
        self.wait()

        distance_line = DashedLine(
            start = bs1[0].get_center(),
            end = ghost.get_center(),
            color = PURPLE
        )
        distance_line_ghost = distance_line.copy().set_opacity(0.6)
        self.play(Create(distance_line))
        self.add(distance_line_ghost)
        self.wait()

        self.play(distance_line.animate.move_to([4.5, -0.3, 0]))
        self.wait()

        distance_brace = BraceLabel(
            distance_line,
            "\\Delta\\vec x = \\vec v \\Delta t",
            brace_direction = DOWN,
        ).set_color(PURPLE)
        self.play(Create(distance_brace))
        self.wait()

class EulerGraphExplanation(MovingCameraScene):
    def construct(self):
        Amp = 2
        dt = 0.1
        k = 5
        m = 1
        period = 2*PI / np.sqrt(k / m)

        ax = Axes(
            x_range = [0, 2*PI],
            y_range = [-2, 2],
            x_length = 7,
            y_length = 3
        )
        labels = ax.get_axis_labels(x_label="t", y_label = "x")

        def spring_func(t):
            return Amp*np.cos(np.sqrt(k/m) * t)
        def derivative(t):
            return -Amp * np.sqrt(k/m) * np.sin(np.sqrt(k/m) * t)
    
        true_func = ax.plot(
            lambda t: spring_func(t),
            color = RED
        ).set_stroke(opacity = 0.5)

        self.play(Create(VGroup(ax, labels)), run_time = 2)
        self.play(Create(true_func), run_time = 1.5)
        self.wait()

        tracker_point = Dot(color = BLUE).move_to(ax.c2p(period+0.2,spring_func(period+0.2),0)).scale(0.9)
        point_coord = MathTex("(t_i, x_i)").scale(0.6).next_to(tracker_point, UR).set_color(BLUE)
        self.play(Create(tracker_point), Create(point_coord))
        self.wait()

        self.play(
            self.camera.frame.animate.move_to(tracker_point.get_center() + DOWN*1).scale(0.5)
        )
        self.wait()

        dt_brace = BraceBetweenPoints(
            ax.c2p(period + 0.2, 0),
            ax.c2p(period + 0.6, 0),
            direction = DOWN
        )
        dt_brace = VGroup(
            dt_brace,
            MathTex("\\Delta t").next_to(dt_brace, DOWN)
        )
        self.play(Create(dt_brace))
        self.wait()

        connecting_line = DashedLine(
            start = tracker_point.get_center(),
            end = ax.c2p(period + 0.2, 0),
            color = PURPLE
        ).set_opacity(0.4)
        self.play(Create(connecting_line))
        self.wait()

        def tangent_func(t):
            return derivative(period + 0.2) * (t - (period + 0.2)) + spring_func(period + 0.2)

        slope_line = ax.plot(
            lambda t: tangent_func(t)
        ).set_color(BLUE_B).set_opacity(0.5)
        self.play(Create(slope_line))
        self.wait()

        self.add(tracker_point.copy())
        self.play(
            tracker_point.animate.move_to(
                ax.c2p(period + 0.6, tangent_func(period + 0.6),0)
            )
        )
        self.wait()

        new_coord = MathTex("(t_{i+1}, x_{i+1})").scale(0.6).next_to(tracker_point, UR).set_color(BLUE)
        new_connecting_line = DashedLine(
            start = tracker_point.get_center(),
            end = ax.c2p(period + 0.6, 0),
            color = PURPLE
        ).set_opacity(0.4)

        self.play(Write(new_coord))
        self.wait()
        self.play(Create(new_connecting_line))
        self.wait()

class AnalyticEulerImplicitSpring(Scene):
    def construct(self):
        Amp = 2
        dt = 0.05
        k = 5
        m = 1

        bs1 = BlockSpring(start = [-3.5,2,0], axis=RIGHT,color=RED, spring_width = 0.5, spring_period=3, length = 5)
        bs2 = BlockSpring(start = [-3.5,0,0], axis=RIGHT,color=BLUE, spring_width = 0.5, spring_period=3, length = 5)
        bs3 = BlockSpring(start = [-3.5,-2,0], axis=RIGHT,color=YELLOW, spring_width = 0.5, spring_period=3, length = 5)
        bs1[0].shift(RIGHT)
        bs2[0].shift(RIGHT)
        bs3[0].shift(RIGHT)

        def pos_func(t):
            return Amp*np.cos(np.sqrt(k/m) * t)

        self.play(Create(VGroup(bs1, bs2, bs3)))
        self.wait()
        time = 4*PI

        #Hooke's Law
        # upward_mover = ValueTracker(bs2[0].get_center()[1])
        # path_head = always_redraw(
        #     lambda: Dot().move_to([bs2[0].get_center()[0], upward_mover.get_value(), 9])
        # )
        # path_tracer = TracedPath(path_head.get_center, dissipating_time=100, stroke_opacity = 0.38, stroke_width=8).set_color(GREEN)
        # connecting_line = always_redraw(
        #     lambda: DashedLine(
        #         start = path_head.get_center(),
        #         end = bs2[0].get_center()
        #     ).set_opacity(0.5)
        # )
        # self.add(connecting_line)
        # self.add(path_head)
        # self.add(path_tracer)

        def euler_step(x, v, dt, k, m):
            a = - (k / m) * x
            v_new = v + a * dt
            x_new = x + v * dt
            return x_new, v_new
        
        def implicit_step(x, v, dt, k, m):
            A = 1 + (dt**2 * k / m)
            B = -dt * k / m

            x_new = (x + dt * v) / A
            v_new = v + dt * (-k / m * x_new)
            
            return x_new, v_new

        x_euler, v_euler = Amp, 0
        x_implicit, v_implicit = Amp, 0
        t = dt
        while t < time + 2*PI:

            x_analytic = pos_func(t)
            x_euler, v_euler = euler_step(x_euler, v_euler, dt, k, m)
            x_implicit, v_implicit = implicit_step(x_implicit, v_implicit, dt, k, m)

            self.play(
                bs1[0].animate.move_to([x_analytic, 2, 0]), 
                bs2[0].animate.move_to([x_euler, 0, 0]), 
                bs3[0].animate.move_to([x_implicit, -2, 0]), 

                rate_func = linear, 
                run_time = dt/2
            )
            t += dt
        self.wait()



class RKGraphExplanation(MovingCameraScene):
    def construct(self):
        Amp = 2
        dt = 0.1
        k = 5
        m = 1
        period = 2*PI / np.sqrt(k / m)

        ax = Axes(
            x_range = [0, 2*PI],
            y_range = [-2, 2],
            x_length = 7,
            y_length = 3
        )
        labels = ax.get_axis_labels(x_label="t", y_label = "x")

        def spring_func(t):
            return Amp*np.cos(np.sqrt(k/m) * t)
        def derivative(t):
            return -Amp * np.sqrt(k/m) * np.sin(np.sqrt(k/m) * t)
    
        true_func = ax.plot(
            lambda t: spring_func(t),
            color = RED
        ).set_stroke(opacity = 0.5)

        self.play(Create(VGroup(ax, labels)), run_time = 2)
        self.play(Create(true_func), run_time = 1.5)
        self.wait()

        tracker_point = Dot(color = BLUE).move_to(ax.c2p(period+0.2,spring_func(period+0.2),0)).scale(0.9)
        point_coord = MathTex("(t_i, x_i)").scale(0.6).next_to(tracker_point, UR).set_color(BLUE)
        self.play(Create(tracker_point), Create(point_coord))
        self.wait()

        self.play(
            self.camera.frame.animate.move_to(tracker_point.get_center() + DOWN*1).scale(0.5)
        )
        self.wait()

        dt_brace = BraceBetweenPoints(
            ax.c2p(period + 0.2, 0),
            ax.c2p(period + 0.6, 0),
            direction = DOWN
        )
        dt_brace = VGroup(
            dt_brace,
            MathTex("\\Delta t").scale(0.7).next_to(dt_brace, DOWN)
        )
        self.play(Create(dt_brace))
        self.wait()

        connecting_line = DashedLine(
            start = tracker_point.get_center(),
            end = ax.c2p(period + 0.2, 0),
            color = PURPLE
        ).set_opacity(0.4)
        self.play(Create(connecting_line))
        self.wait()

        def tangent_func(t):
            return derivative(period + 0.2) * (t - (period + 0.2)) + spring_func(period + 0.2)

        slope_line = ax.plot(
            lambda t: tangent_func(t)
        ).set_color(ORANGE).set_opacity(0.5)
        self.play(Create(slope_line))
        self.wait()
        slope_line_slope = MathTex("m = f(t_i, x_i)").move_to(ax.c2p(1.8, 2.2)).set_color(ORANGE).scale(0.6)
        self.play(Write(slope_line_slope))
        self.wait()

        self.add(ghost := tracker_point.copy().set_opacity(0.4))
        self.play(
            ghost.animate.move_to(
                ax.c2p(period + 0.6, tangent_func(period + 0.6),0)
            )
        )
        self.wait()

        new_coord = MathTex("(t_{i+1}, x_i + k_1)").scale(0.6).next_to(ghost, UR).set_color(BLUE)
        new_connecting_line = DashedLine(
            start = ghost.get_center(),
            end = ax.c2p(period + 0.6, 0),
            color = RED_C
        ).set_opacity(0.4)

        self.play(Write(new_coord))
        self.wait()
        self.play(Create(new_connecting_line))
        self.wait()

        def new_tangent_func(t):
            return derivative(period + 0.6) * (t - (period + 0.6)) + tangent_func(period + 0.6)

        new_slope_line = ax.plot(
            lambda t: new_tangent_func(t)
        ).set_color(GOLD_C).set_opacity(0.5)

        self.play(Create(new_slope_line))
        self.wait()

        new_slope_line_slope = MathTex("m = f(t_{i+1}, x_{i+1})").scale(0.6).next_to(tracker_point, UR).set_color(GOLD_C)
        self.play(FadeTransform(point_coord, new_slope_line_slope))
        self.wait()

        """
        other_func = true_func.copy().set_color(ORANGE).shift(RIGHT*0.15)
        self.play(FadeIn(other_func))
        self.wait()
        self.play(FadeOut(other_func))
        self.wait()

        """

        def between_func(t):
            return (derivative(period + 0.6) + derivative(period + 0.2))/2 * (t - (period + 0.6)) + tangent_func(period + 0.6)
        def between_func_point(t):
            return (derivative(period + 0.6) + derivative(period + 0.2))/2 * (t - (period + 0.2)) + tangent_func(period + 0.2)

        test_line = ax.plot(
            lambda t: between_func(t)
        ).set_color(BLUE_B).set_opacity(0.5)
        self.play(Create(test_line))
        self.wait()

        true_line = ax.plot(
            lambda t: between_func_point(t)
        ).set_color(BLUE_B).set_opacity(0.5)
        self.play(ReplacementTransform(test_line, true_line))
        self.wait()

        self.play(
            tracker_point.animate.move_to(ax.c2p(
                period + 0.6,
                between_func_point(period + 0.6)
            ))
        )

        self.wait()

        newest_coord = MathTex("(t_{i+1}, x_i + (k_1 + k_2)/2").scale(0.5).next_to(tracker_point, RIGHT).set_color(BLUE)
        self.play(Write(newest_coord))
        self.wait()

class RK4GraphExplanationCopy(MovingCameraScene):
    def construct(self):
        ax = Axes(
            x_range = [-.5, 5],
            y_range = [-0.5, 5],
            x_length = 8,
            y_length = 5,
        )

        def func(t):
            return t**2/4

        func_plot = ax.plot(
            lambda t: func(t),
            color = RED
        )
        self.add(VGroup(ax, func_plot))
        self.wait()

        t_labels = VGroup(
            MathTex("t_i").scale(0.6).move_to(ax.c2p(2, -0.5, 0)),
            MathTex("t_i + \\Delta t/2").scale(0.6).move_to(ax.c2p(3, -0.5, 0)),
            MathTex("t_i + \\Delta t").scale(0.6).move_to(ax.c2p(4, -0.5, 0))
        )
        self.play(Write(t_labels))
        self.wait()

        start = ax.c2p(2, func(2))
        middle_1 = ax.c2p(3, 9/4 - 0.5)
        middle_2 = ax.c2p(3, 9/4 + 0.3)
        end = ax.c2p(4, 4.6)

        points = [start, middle_1, middle_2, end]
    
        dots = VGroup(
            *[Dot(color=BLUE).move_to(point) for point in points]
        )

        initial_label = MathTex("(t_i, x_i)").next_to(dots[0], UL).set_color(BLUE)

        connecting_lines = VGroup()
        for i in range(1, len(dots)):
            connecting_lines.add(
                DashedLine(
                    start = dots[0].get_center(),
                    end = dots[i].get_center()
                ).set_opacity(0.5)
            ) 

        line_labels = VGroup(
            MathTex("k_1").scale(0.6).move_to(connecting_lines[0].get_center() + DOWN*0.3),
            MathTex("k_2").scale(0.6).next_to(dots[1], DR, buff = 0.1),
            MathTex("k_3").scale(0.6).next_to(dots[2], UL, buff = 0.1),
            MathTex("k_4").scale(0.6).next_to(dots[3], UL, buff = 0.1)
        )

        slope_lines = VGroup(connecting_lines[0].copy())

        for i in range(2):
            slope_lines.add(
                Line(
                    start = dots[0].get_center(),
                    end = dots[i+2]
                ).scale(0.38**i/2).move_to(dots[i+1].get_center())
            )
        slope_lines.add(
            Line(
                start = dots[0].get_center(),
                end = dots[1].get_center() + DOWN*0.1
            ).scale(0.6).move_to(dots[3].get_center())
        )

        final_point = Dot(color = BLUE).scale(1.2).move_to(ax.c2p(4, func(4) - 0.09))
        final_connection = DashedLine(
            start = dots[0].get_center(),
            end = final_point.get_center()
        ).set_opacity(0.7)
        final_label = MathTex("(t_{i+1}, x_{i+1})").next_to(final_point, DR).set_color(BLUE)


        #Animating
        self.play(Create(dots[0]), Write(initial_label))
        self.wait()
        self.play(Create(connecting_lines[0]), Write(line_labels[0]))

        for i in range(3):
            self.play(Create(dots[i+1]))
            self.wait()

            self.play(Create(slope_lines[i+1]), Write(line_labels[i+1]))
            self.wait()

            if i < 2:
                self.play(
                    Transform(slope_lines[i+1].copy(), connecting_lines[i+1])
                )
                self.wait()
            else:
                continue
            

        self.play(
            *[Transform(mob, final_connection) for mob in connecting_lines]
        )
        self.wait()
        self.play(Create(final_point), Write(final_label))
        self.wait()

        final_formula = VGroup(
            Text("RK4"),
            MathTex("k_1 = f(t_i, x_i)"),
            MathTex("k_2 = f(t_i + \\Delta t/2, x_i + k_1\\Delta t/2)"),
            MathTex("k_3 = f(t_i + \\Delta t/2, x_i + k_2\\Delta t/2)"),
            MathTex("k_4 = f(t_i + \\Delta t, x_i + k_3\\Delta t)"),
            MathTex("x_{i+1} = \\frac{\\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)")
        ).arrange(DOWN).scale(0.58).move_to(ax.c2p(1.7, 4.5))
        final_formula.add(SurroundingRectangle(final_formula, color = GOLD, buff = 0.05))

        self.play(Write(final_formula))
        self.wait()
