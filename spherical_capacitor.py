from manim import *
import numpy as np

class FinalAnimation(Scene):
    def construct(self):
        None

#FADE OUT EVERYTHING ON SCREEN
def fade_all(self):
    self.play(
        *[FadeOut(mob)for mob in self.mobjects]
    )

#Spherical Capacitor Diagram
inner_radius = ValueTracker(2.3)
outer_radius =  ValueTracker(3.3)

inner_shell = Circle(
    radius = inner_radius.get_value(),
    color = RED
)
inner_theta = PI/6
inner_radius_vector = Vector(color=RED)
inner_radius_vector.put_start_and_end_on(
    inner_shell.get_center(),
    RIGHT*inner_radius.get_value()*np.cos(inner_theta) + UP*inner_radius.get_value()*np.sin(inner_theta)
)
a = MathTex("a").next_to(inner_radius_vector.get_center(), UP, buff=0.1).scale(0.8)
a.set_color(inner_radius_vector.get_color())
Qpos = MathTex("Q").next_to(inner_shell, UP, buff = -0.6).set_color(inner_shell.get_color())

outer_shell = Circle(
    radius = outer_radius.get_value(),
    color = PURPLE
)
outer_theta = PI/3
outer_radius_vector = Vector(color=PURPLE)
outer_radius_vector.put_start_and_end_on(
    outer_shell.get_center(),
    RIGHT*outer_radius.get_value()*np.cos(outer_theta) + UP*outer_radius.get_value()*np.sin(outer_theta)
)
b = MathTex("b").next_to(outer_radius_vector.get_center(), UP, buff=0.2).scale(0.8)
b.set_color(outer_radius_vector.get_color())
Qneg = MathTex("-Q").next_to(outer_shell, UP, buff = -0.6).set_color(outer_shell.get_color())

capacitor = VGroup(a, b, inner_shell, outer_shell, inner_radius_vector, outer_radius_vector, Qpos, Qneg)

class Title(Scene):
    def construct(self):
        self.wait(2)
        question = Text("""
                        Find the capacitance of two
                        concentric spherical shells. 
                        The inner shell has radius a
                        and charge +Q, the outer shell
                        has radius b > a and charge -Q.""")
        
        self.play(Create(question), run_time = 10)
        self.wait(3)


        self.play(FadeOut(question), run_time=0.3)

        self.play(Create(capacitor))
        self.play(
            capacitor.animate.scale(0.8).shift(RIGHT*3.6),
            FadeIn(question.shift(LEFT*2.9).scale(0.8)),
            run_time=1.5
        )

        self.play(question.animate.shift(DOWN*1))
        self.wait(3)

        #Credit to FillahAlamsyah from https://github.com/FillahAlamsyah/Manim-MathAnimation/blob/7543445378ad9ba5a5730955724de5914c28b602/Countdown.py
        number = Text("5").scale(2).next_to(question.get_center(), UP, buff=2.5)
        self.play(Write(number))
        for i in range(4,-2,-1):
            if i != -1 :
                self.play(Transform(number,
                Text(str(i)).next_to(question.get_center(), UP, buff=2.5).scale(2)),run_time=1)
                self.wait(0.7)
            elif i == -1 :
                self.play(Transform(number,
                Text(str(i+1)).next_to(question.get_center(), UP, buff=2.5).scale(2)),run_time=0.05)
                self.play(FadeOut(number), run_time=0.5)


        self.wait(5)
        self.play(FadeOut(question))

        #PROBLEM SOLVING
        self.wait(5)
        final_equation = MathTex("C = \\frac{Q}{V}").shift(LEFT*3+UP*2.5).scale(1.3)
        V = MathTex("V")
        potential_definition = MathTex("V = -\\int_b^a \\vec{E}","\\cdot \\vec{ds}").shift(final_equation.get_center()+DOWN*2)
        potential_case_specific = MathTex("V = -\\int_b^aE","dr").shift(final_equation.get_center()+LEFT*0.2+DOWN*2)
        E = MathTex("\\vec{E}")
        electric_field = MathTex("")

        gaussian_sphere = DashedVMobject(Circle(
            radius = ((inner_radius.get_value() + outer_radius.get_value())/2)*0.8,
            color = WHITE,
        ).move_to(inner_shell.get_center())
        )
        gaussian_sphere.set_fill(opacity=0.7)
        

        gauss_law_definition = MathTex("\\oint_S \\vec{E} \\cdot d\\vec{A} = \\frac{Q_{enc}}{\\epsilon_0}")
        gauss_law_specific = MathTex("\\iint_S E dA = \\frac{Q_{enc}}{\\epsilon_0}")
        self.play(Create(final_equation))
        self.wait(6)
        self.play(Create(potential_definition))
        self.wait(10)
        self.play(FadeTransform(potential_definition, potential_case_specific))
        self.wait(6)
        self.play(Create(E.shift(final_equation.get_center()+DOWN*3.5)))
        self.wait(6)
        self.play(Create(gauss_law_definition.shift(final_equation.get_center()+LEFT*0.2+DOWN*5)))
        self.wait(6)
        self.play(Create(gaussian_sphere))
        capacitor.add(gaussian_sphere)
        self.wait(6)
        self.play(FadeTransformPieces(gauss_law_definition,gauss_law_specific.shift(final_equation.get_center()+LEFT*0.2+DOWN*5)))
        self.wait()

        self.play(
            final_equation.animate.shift(LEFT*2).scale(0.9),
            potential_case_specific.animate.shift(UP*2+RIGHT*1.9),
            gauss_law_specific.animate.shift(RIGHT*6 + UP*4.9),
            FadeOut(E),
            capacitor.animate.shift(DOWN*1).scale(0.9)
        )
        self.wait(4)

        rect = SurroundingRectangle(gauss_law_specific)
        self.play(Create(rect))
        self.play(Uncreate(rect.reverse_direction()))
        rect = SurroundingRectangle(potential_case_specific)
        self.play(Create(rect))
        self.play(Uncreate(rect.reverse_direction()))
        self.wait(1)
        rect = SurroundingRectangle(final_equation)
        self.play(Create(rect))
        self.play(Uncreate(rect.reverse_direction()))
        self.wait(1)

        retention_group = VGroup(final_equation, capacitor, gauss_law_specific, potential_case_specific)
        #TRANSITION TO FINDING E FIELD
        temp_gauss_specific = gauss_law_specific
        everything_except_gauss_law_specific = Group()
        for mob in self.mobjects:
            if mob != temp_gauss_specific:
                everything_except_gauss_law_specific.add(mob)
        self.play(
            temp_gauss_specific.animate.move_to([0,0,0]),
            *[FadeOut(mob)for mob in everything_except_gauss_law_specific],
        )

        self.wait(2)
        self.remove(temp_gauss_specific)
        self.wait(2)
        self.play(
            *[FadeIn(mob) for mob in retention_group],
            gauss_law_specific.animate.move_to([3.5,2.5,0])
        )
        self.wait(3)

        #TRANSITION TO FINDING POTENTIAL
        temp_potential_specific = potential_case_specific
        everything_except_potential_case_specific = Group()
        for mob in self.mobjects:
            if mob != temp_potential_specific:
                everything_except_potential_case_specific.add(mob)
        self.play(
            temp_potential_specific.animate.move_to([0,0,0]),
            *[FadeOut(mob)for mob in everything_except_potential_case_specific],
        )
        self.wait(2)
        self.remove(temp_potential_specific)
        self.wait(2)
        self.play(
            *[FadeIn(mob) for mob in retention_group],
            potential_case_specific.animate.move_to([-1,2.5,0])
        )
        self.wait(3)

        final_potential = MathTex("V = \\frac{Q}{4\\pi\\epsilon_0}\\frac{b-a}{ab}").move_to([-1,2.5,0])
        self.play(ReplacementTransform(potential_case_specific, final_potential))

        capacitance_equation = MathTex("C = 4\\pi\\epsilon_0 \\left(\\frac{ab}{b-a}\\right)")
        self.wait(2)
        self.play(Write(capacitance_equation.shift(LEFT*2.3+DOWN*0.5).scale(1.5)))
        self.play(Create(rect))
        self.wait(10)
        fade_all(self)
        

class FindingEField(Scene):
    def construct(self):
        starting_equation = MathTex(
            "\\iint_S",          #0
            "E",                 #1
            "dA",                #2
            "=",                 #3
            "{Q_{enc}",          #4
            "\\over",            #5
            "\\epsilon_0}"       #6
        )
        gaussian_sphere = DashedVMobject(Circle(
            radius = ((inner_radius.get_value() + outer_radius.get_value())/2),
            color = WHITE,
        ).move_to(inner_shell.get_center())
        )
        gaussian_sphere.set_fill(opacity=0.7)
        capacitor.add(gaussian_sphere)
        self.add(starting_equation)
        self.play(FadeIn(capacitor.shift(RIGHT*4.5, UP*1.5).scale(0.6)))
        self.wait(5)
        self.play(Swap(starting_equation[0], starting_equation[1]))
        self.wait(4)

        surface_area = MathTex("4\\pi r^2")
        self.play(
            ReplacementTransform(
                Group(starting_equation[0], starting_equation[2]),
                surface_area.shift(LEFT*0.3+UP*0.05)
            ),
            starting_equation[1].animate.shift(RIGHT*0.4)
        )

        self.wait(4)
        Qover = MathTex("Q\\over").shift(RIGHT*1.2+UP*0.25)
        self.play(
            ReplacementTransform(starting_equation[4:6], Qover),
            starting_equation[3].animate.shift(RIGHT*0.1),
            starting_equation[6].animate.shift(LEFT*0.1)
        )
        self.wait(3)
        substepgroup1 = VGroup(starting_equation[1], starting_equation[3], surface_area, starting_equation[6], Qover)
        self.play(substepgroup1.animate.shift(DOWN*2))

        E_field = MathTex("\\vec{E} = \\frac{1}{4\\pi\\epsilon_0}\\frac{Q}{r^2}\\hat{r}")
        self.play(Write(E_field.scale(1.2)))
        self.wait(10)
        fade_all(self)


class FieldtoPotential(Scene):
    def construct(self):

        E_field = MathTex(
            "\\vec{E}",     #0
            "=",            #1
            "{1",           #2
            "\\over",       #3
            "4",            #4
            "\\pi",         #5
            "\\epsilon_0}", #6
            "{Q",           #7
            "\\over",       #8
            "r",            #9
            "^2}",          #10
            "\\hat{r}"      #11
        )
        gaussian_sphere = DashedVMobject(Circle(
            radius = ((inner_radius.get_value() + outer_radius.get_value())/2),
            color = WHITE,
        ).move_to(inner_shell.get_center())
        )
        capacitor.add(gaussian_sphere)
        
        potential = MathTex(
            "V",            #0
            "=",            #1
            "-",            #2
            "\\int_b^a",    #3
            "E",            #4
            "dr"            #5
        )

        self.add(potential)

        self.play(
            FadeIn(capacitor.shift(RIGHT*4.5, UP*1.5).scale(0.6)),
            run_time=0.5
        )
        self.wait(3)

        self.play(
            Write(E_field),
            E_field.animate.shift(DOWN*2)
        )

        self.wait(4)
        E = MathTex("E")
        self.play(
            FadeTransform(E_field[0], E.shift(DOWN*2+LEFT*1.25)),
            FadeOut(E_field[11])
        )
        self.wait(2)
        self.play(
            FadeOut(potential[4]),
            FadeOut(E_field[1]),
            FadeOut(E),
            E_field[2:11].animate.shift(UP*2+RIGHT*0.5),
            potential[0:4].animate.shift(LEFT*0.5),
            potential[5].animate.shift(RIGHT*0.75)
        )

        final_potential_integral = VGroup(
            E_field[2:11],
            potential[0:4],
            potential[5]
        )

        self.wait()
        self.play(final_potential_integral.animate.shift(UP*2.5+LEFT*2))
        self.wait()
        self.play(
            Swap(E_field[2:7], potential[3]),
            Group(
                E_field[7:11],
                potential[0:3],
                potential[5]
            ).animate.shift(0.1*LEFT+UP*0.05),   
        )
        self.wait(3)
        self.play(Swap(E_field[7], E_field[2]))
        r_to_minus_2 = MathTex("r^{-2}").shift(UP*2.6+LEFT*0.65)
        self.wait(2)
        self.play(ReplacementTransform(VGroup(E_field[2],E_field[8:11]), r_to_minus_2))
        self.wait(2)
        integral_evaluated = MathTex(
            "=",            #0
            "-",            #1
            "{Q",           #2
            "\\over",       #3
            "4",            #4
            "\\pi",         #5
            "\\epsilon_0}", #6
            "\\left(",      #7
            "-",            #8
            "r",            #9
            "^{-1}",        #10
            "]"             #11    
            "^a",           #12
            "_b",           #13
            "\\right)"      #14
        ).shift(LEFT*1.65+UP*1)
        self.play(Write(integral_evaluated))
        self.wait()
        self.play(
            FadeOut(integral_evaluated[1]),
            FadeOut(integral_evaluated[8]),
            integral_evaluated[2:8].animate.shift(LEFT*0.4),
            integral_evaluated[9:15].animate.shift(LEFT*.8)
        )
        self.wait(3)

        values_plugged_in = MathTex(
            "\\left(",   #0
            "{1",       #1
            "\\over",   #2
            "a}",       #3
            "-",        #4
            "{1",       #5
            "\\over",   #6
            "b}",       #7
            "\\right)"  #8
        ).shift(LEFT*1.05+UP*1)

        self.play(ReplacementTransform(Group(integral_evaluated[7],integral_evaluated[9:14]), values_plugged_in))
        algebraically_rearranged = MathTex(
            "V = ",
            "\\frac{Q}{4\\pi\\epsilon_0}",
        "\\left(",          #0
            "{b",           #1
            "-",            #2
            "a",            #3
            "\\over",       #4
            "ab}",          #5
            "\\right)"      #6
        ).shift(LEFT*2.1+DOWN*0.5)
        self.wait(3)
        self.play(Write(algebraically_rearranged))
        self.wait(3)

        self.play(
            *[mob.animate.shift(LEFT*2.2) for mob in self.mobjects],
        )
        self.play(
            capacitor.animate.shift(DOWN*1.5),
        )
        self.play(capacitor.animate.scale(1.5))
        inner_shell_group = VGroup(inner_shell, inner_radius_vector, a, Qpos)
        outer_shell_group = VGroup(outer_shell, outer_radius_vector, b, Qneg)

        self.wait(2)
        self.play(inner_shell_group.animate.scale(0.8))
        self.wait()
        self.play(outer_shell_group.animate.scale(1.1))
        self.wait(3)
        self.play(
            inner_shell_group.animate.scale(1.8),
            outer_shell_group.animate.scale(0.7)
        )
        self.play(
            inner_shell_group.animate.scale(0.6),
            outer_shell_group.animate.scale(1.5)
        )

        self.wait(4)
        self.play(
            capacitor.animate.scale(0.8)
        )

        self.wait(4)
        fade_all(self)

class SphereAnimation(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=60*DEGREES, theta=-45*DEGREES)
        inner_shell3D = Sphere(
            radius = 2.3
        ).set_color(RED)
        inner_shell3D.set_opacity(0.4)
        outer_shell3D = Sphere(
            radius = 3.3,
        ).set_color(PURPLE)
        outer_shell3D.set_opacity(0.3)


        self.begin_ambient_camera_rotation(
            rate = PI/10, about="theta" #per second
        )
        self.play(Create(inner_shell3D))
        self.play(Create(outer_shell3D))

        self.wait(6)
