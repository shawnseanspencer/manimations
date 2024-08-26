from manim import *

class AchillesTortoise(MovingCameraScene):
    def construct(self):
        tortoise = ImageMobject("/Users/iansackin/Manim/Pics/Tortoise.jpeg").shift(RIGHT*4.5)
        achilles = ImageMobject("/Users/iansackin/Manim/Pics/Achilles.png").shift(LEFT*4.5).scale(0.5).set_y(tortoise.get_y())

        connecting_line = Line(
            end = tortoise.get_center(),
            start = achilles.get_center() + RIGHT*0.2
        ).set_opacity(0.6)

        self.play(FadeIn(Group(tortoise, achilles, connecting_line)))
        self.wait()

        notches = VGroup()
        length = connecting_line.get_length()
        for i in range(1, 10):
            notches.add(
                Line(
                    start = [0,0,0],
                    end = [0,0.6,0],
                    stroke_width = 1
                ).move_to(connecting_line.get_end() + LEFT*length/2**i).set_opacity(0.8)
            )

        self.play(Create(notches))
        self.wait()

        for i in range(len(notches) - 3):
            self.play(achilles.animate.set_x(notches[i].get_x()).shift(0.22*LEFT))
            self.wait()

        self.play(
            tortoise.animate.scale(0.05),
            achilles.animate.scale(0.05).set_x(notches[-3].get_x()).shift(0.25*0.05*LEFT),
            self.camera.frame.animate.scale(0.05, about_point = connecting_line.get_end())
        )
        self.wait()

class MegumiBum(MovingCameraScene):
    def construct(self):
        ax = Axes(
            x_range = [0, 10],
            y_range = [0, 10]
        )

        val = 5.3

        def f(x):
            return 0.9 * x
    

        func = ax.plot(lambda x: f(x)).set_color(BLUE)
        labels = ax.get_axis_labels(x_label="\\text{Megumi}", y_label = "\\text{Bum}")

        self.play(Create(VGroup(ax, func, labels)))
        self.wait()

        me = Dot(color = RED).move_to(
            ax.c2p(val, 0)
        ).scale(1.5)

        px = DashedLine(
            stroke_width = 3,
            start = me.get_center(),
            end = ax.c2p(val, f(val))
        ).set_opacity(0.3)
        py = DashedLine(
            stroke_width = 3,
            start = ax.c2p(0, f(val)),
            end = ax.c2p(val, f(val))
        ).set_opacity(0.3)

        self.play(Create(me))
        self.wait()
        self.play(Create(px))
        self.wait()
        self.play(Create(py), reversed = True)
        self.wait()

        self.play(me.animate.move_to(px.get_end()))

        self.wait()

        epsilon = ValueTracker(0.7)
        delta = ValueTracker(0.6)

        d_lines = VGroup(
            always_redraw(
                lambda: DashedLine(
                    color = GREEN,
                    start = ax.c2p(val + delta.get_value(),0),
                    end = ax.c2p(val + delta.get_value(), 10)
                ).set_opacity(0.7)
            ),
            always_redraw(
                lambda: DashedLine(
                    color = GREEN,
                    start = ax.c2p(val - delta.get_value(),0),
                    end = ax.c2p(val - delta.get_value(), 10)
                ).set_opacity(0.7)
            )
        )

        e_lines = VGroup(
            always_redraw(
                lambda: DashedLine(
                    color = PURPLE,
                    start = ax.c2p(0, f(val) + epsilon.get_value()),
                    end = ax.c2p(10, f(val) + epsilon.get_value())
                ).set_opacity(0.8)
            ),
            always_redraw(
                lambda: DashedLine(
                    color = PURPLE,
                    start = ax.c2p(0, f(val) - epsilon.get_value()),
                    end = ax.c2p(10, f(val) - epsilon.get_value())
                ).set_opacity(0.8)
            )
        )

        self.play(Create(d_lines))
        self.wait()
        self.play(delta.animate.set_value(0.3))
        self.wait()

        self.play(Create(e_lines))
        self.wait()
        self.play(epsilon.animate.set_value(0.4))

        self.wait()

        self.play(
            epsilon.animate.set_value(1),
            delta.animate.set_value(1)
        )
        self.wait()

        coord = MathTex("(53.245\\%, b)").scale(0.7).next_to(me, UL, buff = 0.2).set_color(RED)
        self.play(Write(coord))
        self.wait()

        for _ in range(4):
            temp_epsilon = epsilon.get_value()
            temp_delta = delta.get_value()

            self.play(
                epsilon.animate.set_value(temp_epsilon/2),
                delta.animate.set_value(temp_delta/2)
            )

        self.wait()

        self.play(
            epsilon.animate.set_value(1.5),
            delta.animate.set_value(1.5)
        )

        e_brace = always_redraw(lambda:  BraceBetweenPoints(
            ax.c2p(3.5, f(val)),
            ax.c2p(3.5, f(val) + epsilon.get_value()),
            direction = LEFT,
            color = PURPLE
        ))
        d_brace = always_redraw(lambda:  BraceBetweenPoints(
            ax.c2p(val, 2.5),
            ax.c2p(val + delta.get_value(), 2.5),
            direction = DOWN,
            color = GREEN
        ))
        e_brace_group = VGroup(
            e_brace, always_redraw(lambda: MathTex("\\epsilon").move_to(e_brace.get_center() + LEFT*0.5).set_color(PURPLE))
        )
        d_brace_group = VGroup(
            d_brace, always_redraw(lambda: MathTex("\\delta").move_to(d_brace.get_center() + DOWN*0.5).set_color(GREEN))
        )

        self.play(Write(e_brace_group))
        self.wait()
        self.play(Write(d_brace_group))
        self.wait()

        self.wait()
        hole = Circle(radius = 0.1).set_fill(BLACK, opacity = 1).set_stroke(RED).move_to(me.get_center())
        self.play(me.animate.shift(UP*1.5))
        self.wait()
        self.play(FadeIn(hole))
        self.wait()

        self.play(
            epsilon.animate.set_value(0.15),
            delta.animate.set_value(0.15),
            run_time = 5,
            rate_fun = linear
        )
        self.wait()

class AlgebraStuff(Scene):
    def construct(self):

        sum = VGroup(
            MathTex("\\frac{1}{2} + \\frac{1}{4} + \\frac{1}{8} + \\frac{1}{16} + ..."),
            MathTex("\\lim_{n\\to \\infty}","\\sum_{i=1}^{n} \\frac{1}{2^i}", "=", "\\frac{1}{2} + \\frac{1}{2^2} + \\frac{1}{2^3} + ...")
        ).arrange(DOWN)
        self.play(Write(sum[0]))
        self.wait()

        self.play(Write(sum[1][1]))
        self.wait()
        self.play(Write(sum[1][2:]))
        self.wait()

        rects = VGroup(
            SurroundingRectangle(sum[0], buff = 0.1),
            SurroundingRectangle(sum[1][3], buff = 0.2)
        )

        self.play(Create(rects))
        self.wait()
        self.play(Write(sum[1][0]))
        self.wait()

class AchillesReachesTortoise(Scene):
    def construct(self):

        t_to_infty = MathTex("n \\to \\infty").shift(UP*1.5).scale(1.3)

        tortoise = ImageMobject("/Users/iansackin/Manim/Pics/Tortoise.jpeg").shift(RIGHT*4.5)
        achilles = ImageMobject("/Users/iansackin/Manim/Pics/Achilles.png").shift(LEFT*4.5).scale(0.5).set_y(tortoise.get_y())

        connecting_line = Line(
            end = tortoise.get_center(),
            start = achilles.get_center() + RIGHT*0.2
        ).set_opacity(0.6)

        self.add(Group(tortoise, achilles, connecting_line, t_to_infty))

        self.play(achilles.animate.shift(RIGHT*connecting_line.get_length()))
        self.wait()

class Statements(Scene):
    def construct(self):
        limit = VGroup(
            Tex("The limit as $x$ approaches 53.245 of Bum$(x)$ is b"),
            MathTex("\\lim_{x\\to 53.245} \\text{Bum}(x) = b")
        ).arrange(DOWN)
        self.play(Write(limit[0]))
        self.wait()
        self.play(Write(limit[1]))
        self.wait()

        self.play(FadeOut(*self.mobjects))
        self.wait()

        statement = VGroup(
            Tex("All we care about is what value Bum$(x)$"),
            Tex("approaches as $x$ gets arbitrarily close to 53.245")
        ).arrange(DOWN)
        self.play(Create(statement), run_time = 4)
        self.wait()

class DistanceSubdivision(Scene):
    def construct(self):
        distance = Line(
            start = [-4, 0, 0],
            end = [4, 0, 0]
        )

        self.add(distance)

        first = VGroup()

        for i in range(6):
            for j in range(1, 4**i):
                dx = distance.get_length()/(4**i)
                first.add(
                    Line(
                        start = [0, 0, 0],
                        end = [0,0.8,0],
                        stroke_width = 0.8,
                    ).move_to(distance.get_start() + j * dx*RIGHT) 
                )

            self.play(Create(first.copy()))

        self.wait()


class EpsilonDelta(Scene):
    def construct(self):
        epsilon_constraint = MathTex("|\\text{Bum}(x)} - b| < \\epsilon").set_color(RED)

        self.add(epsilon_constraint)

class InfinityExplanation(MovingCameraScene):
    def construct(self):
        ax = Axes(
            x_range = [0, 10],
            y_range = [0, 10]
        )

        def f(x):
            return 0.9 * x
    

        func = ax.plot(lambda x: f(x)).set_color(BLUE)
        temp = ax.get_axis_labels(x_label="\\text{Megumi}", y_label = "\\text{Bum}")

        self.play(Create(VGroup(ax, func, temp)))
        self.wait()

        labels = ax.get_axis_labels(x_label="t", y_label = "s")
        self.play(ReplacementTransform(temp[0], labels[0]), ReplacementTransform(temp[1], labels[1]))
        self.wait()

        gojo_line = ax.plot(
            lambda x: 7.3,
            x_range = [0, 15.667],
            color = WHITE
        )
        gojo_label = ImageMobject("/Users/iansackin/Manim/Pics/nah_id_win.jpg").scale(0.3).next_to(gojo_line, UP).shift(LEFT*4)

        infty_line = ax.plot(
            lambda x: 5.1,
            x_range = [0, 15.667],
            color = LIGHT_GRAY
        ).set_opacity(0.7)
        infty_label = Text("Infinity Barrier").scale(0.7).next_to(infty_line, UP).shift(LEFT*4).set_color(LIGHT_GRAY)

        tracker_point = Dot(color = RED).move_to(ax.c2p(0,0))
        self.play(tracker_point.animate.move_to(ax.c2p(10,f(10))), rate_func = linear, run_time = 3)
        tracker_point.move_to(ax.c2p(0,0))
        self.wait()

        self.play(FadeIn(Group(gojo_line, gojo_label)))
        self.wait()
        self.play(Create(VGroup(infty_line, infty_label)))
        self.wait()


        asymptote = VGroup(
            ax.plot(lambda x: f(x), x_range = [0, 5.667]).set_color(BLUE),
            ax.copy().shift(
                ax.c2p(5.667, 0) - ax.c2p(0, 0)
            ).plot(lambda x: -1/(1.55115**(x-5.459)*5) + 7.3).set_color(BLUE)
        )

        self.play(ReplacementTransform(func, asymptote))
        self.wait()

        self.play(self.camera.frame.animate.shift(RIGHT*5), run_time = 4)
        self.wait()

        self.play(tracker_point.animate.move_to(ax.c2p(5.667,f(5.667))), rate_func = linear, run_time = 1.5)
        self.play(MoveAlongPath(tracker_point, asymptote[1]), rate_func = linear, run_time = 3)
        self.wait()


class TimeSpaceConstant(ThreeDScene):
    def construct(self):

        self.set_camera_orientation(
            zoom = 1.3
        )

        length = 10
        start = [-length/2, 0, 0]

        lines = VGroup()
        for i in range(1, 15):
            lines.add(
                Line(
                    start = start,
                    end = start + 1/2**i * length * RIGHT
                )
            )

            start += 1/2**i * length * RIGHT

        self.play(Create(lines))
        self.wait()

        parts = VGroup(*lines.copy()).arrange(DOWN).shift(DOWN * 1)
        self.play(ReplacementTransform(lines, parts))
        self.wait()

        statement = Text("These all take the same time to traverse").next_to(parts, UP).scale(0.9)
        self.play(Write(statement))
        self.wait()
