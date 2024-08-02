from manim import *
import numpy as np
from math import *

#Sets the color palatte for the video. Can access through dictionary
#or just use as a reference. I do both regularly
PALATTE = {
    "Titles": BLUE_C,
    "Content": BLUE_D,
    "Gray": GRAY
}

def fade_all(self):
    self.play(
        *[FadeOut(mob) for mob in self.mobjects]
    )


#Flashes mobject a color for an alloted time, very common animation
#so it was useful to just write a function
def highlight(self, mob, color, rt=1) -> None:
    rect = SurroundingRectangle(mob, color)
    self.play(FadeIn(rect), run_time = rt)
    return rect

def get_1D_lattice(
    point_count_h, 
    shape = Cube,
    side_length = 0.2,
    i_ah = 1,
    real_color = BLUE
):
    lattice_1D = VGroup()
    ah_lattice = ValueTracker(i_ah) #Lattice constant horizontal

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

    return lattice_1D
    

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
    lattice_2D.move_to([0,0,0])

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

class DiffractionPatterns(ThreeDScene):
    def construct(self):


        """
        vector_conditions = VGroup()
        vector_conditions.add(Tex(
            "$\\mathbf{a}$",    #0
            "$_i$",             #1
            "$\\mathbf{b}$",    #2
            "$_j$",             #3
            "$=$",              #4
            "$2\\pi n$",        #5
            "$\\delta$",        #6
            "$_i$",             #7
            "$_j$"              #8
        ))
        VGroup(vector_conditions[0][0:2], vector_conditions[0][7]).set_color(BLUE)
        VGroup(vector_conditions[0][2:4], vector_conditions[0][8]).set_color(RED)
        vector_conditions.add(SurroundingRectangle(vector_conditions, PURPLE_E, buff = 0.2))
        vector_conditions.add(BackgroundRectangle(vector_conditions, buff = 0.2))
        vector_conditions = VGroup(vector_conditions[2], vector_conditions[1], vector_conditions[0])
        self.play(Write(vector_conditions))
        self.wait()
        self.play(vector_conditions.animate.shift(UP*3.2 ))
        self.add_fixed_in_frame_mobjects(vector_conditions)
        self.wait()
        """
        self.set_camera_orientation(gamma = 90*DEGREES)

        lattice_1D = VGroup()
        ah_lattice = ValueTracker(1) #Lattice constant horizontal
        av_lattice = ValueTracker(1) #Lattice constant vertical

        lattice_center = [
            ValueTracker(0),
            ValueTracker(0),
            ValueTracker(0),
        ]
        point_count_h = 5 #Horizontal point count
        for i in range(point_count_h):
            start = -point_count_h/2 * ah_lattice.get_value()
            point = Cube(
                side_length = 0.2
            ).set_color(BLUE).move_to(RIGHT*(start + ah_lattice.get_value()*i))
            lattice_1D.add(point)
        lattice_1D = lattice_1D.move_to([0,0,0])

        lattice_2D = VGroup()
        point_count_v = 5

        for i in range(point_count_v):
            start = -point_count_v/2 * av_lattice.get_value()
            lattice_2D.add(lattice_1D.copy().move_to(DOWN*(start + av_lattice.get_value()*i)))
        lattice_2D.move_to([0,0,0])

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
            point = Cube(
                side_length = 0.2
            ).set_color(RED).move_to(RIGHT*(start + 1/ah_lattice.get_value()*i))
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

        
        self.play(Write(lattice_2D))
        self.wait()
        self.play(Write(recip_2D))
        self.wait()
        real_label = always_redraw(
            lambda: Tex("Real Space").next_to(lattice_2D, RIGHT, buff=-0.5).set_color(BLUE).rotate(angle=-90*DEGREES)
        ) #Real Space Label
        recip_label = always_redraw(
            lambda: Tex("Reciprocal Space").next_to(recip_2D, RIGHT, buff=-1.1).set_color(RED).rotate(angle=-90*DEGREES)
        ) #Reciprocal Space Label
        lattice_bracket_h = always_redraw(
            lambda: BraceBetweenPoints(lattice_2D[0][1].get_center(), lattice_2D[0][2].get_center()).set_color(GREEN)
        )
        recip_bracket_h = always_redraw(
            lambda: BraceBetweenPoints(recip_2D[0][1].get_center(), recip_2D[0][2].get_center()).set_color(LOGO_GREEN)
        )
        lattice_bracket_v = always_redraw(
            lambda: BraceBetweenPoints(lattice_2D[0][2].get_center(), lattice_2D[1][2].get_center()).set_color(PURPLE)
        )
        recip_bracket_v = always_redraw(
            lambda: BraceBetweenPoints(recip_2D[0][2].get_center(), recip_2D[1][2].get_center()).set_color(PURPLE_E)
        )
        """
        ah_counter = always_redraw(
            lambda: DecimalNumber(ah_lattice.get_value()).set_color(GREEN).move_to(lattice_bracket_h.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        av_counter = always_redraw(
            lambda: DecimalNumber(av_lattice.get_value()).set_color(PURPLE).move_to(lattice_bracket_v.get_tip()+RIGHT*0.2).rotate(angle=-90*DEGREES, axis=OUT)
        )
        recip_ah_counter = always_redraw(
            lambda: DecimalNumber(1/ah_lattice.get_value()).set_color(LOGO_GREEN).move_to(recip_bracket_h.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        recip_av_counter = always_redraw(
            lambda: DecimalNumber(1/av_lattice.get_value()).set_color(PURPLE_E).move_to(recip_bracket_v.get_tip()+RIGHT*0.2).rotate(angle=-90*DEGREES, axis=OUT)
        )
        """

        ah_counter = always_redraw(
            lambda: MathTex("L_x").set_color(GREEN).move_to(lattice_bracket_h.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        av_counter = always_redraw(
            lambda: MathTex("L_y").set_color(PURPLE).move_to(lattice_bracket_v.get_tip()+RIGHT*0.2).rotate(angle=-90*DEGREES, axis=OUT)
        )
        recip_ah_counter = always_redraw(
            lambda: MathTex("2\\pi / L_x").set_color(LOGO_GREEN).move_to(recip_bracket_h.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        recip_av_counter = always_redraw(
            lambda: MathTex("2\\pi / L_y").set_color(PURPLE_E).move_to(recip_bracket_v.get_tip()+RIGHT*0.2).rotate(angle=-90*DEGREES, axis=OUT)
        )


        brackets = VGroup(lattice_bracket_h, lattice_bracket_v, recip_bracket_h, recip_bracket_v)
        twoD_counters = VGroup(ah_counter, av_counter, recip_ah_counter, recip_av_counter)
        self.play(Write(VGroup(real_label, recip_label)))
        self.wait()
        self.play(Write(brackets))
        self.wait()
        self.play(Write(twoD_counters))
        self.wait()

        self.play(
            ah_lattice.animate.set_value(0.8),
            av_lattice.animate.set_value(0.8),
            lattice_center[1].animate.set_value(-3),
            recip_center[1].animate.set_value(3),
        )
        self.wait()
        self.play(Circumscribe(recip_bracket_h, color=GOLD), Circumscribe(lattice_bracket_h, color=GOLD), run_time = 2)
        self.wait()
        self.play(Circumscribe(recip_bracket_v, color=GOLD), Circumscribe(lattice_bracket_v, color=GOLD), run_time = 2)
        self.wait()
        self.play(
            ah_lattice.animate.set_value(1.2),
            av_lattice.animate.set_value(1.2),
        )
        self.wait()
        self.play(ah_lattice.animate.set_value(0.7))

        
        self.move_camera(phi = -50*DEGREES)
        self.play(
            recip_center[2].animate.set_value(-4),
            lattice_center[2].animate.set_value(4),
            recip_center[1].animate.set_value(0),
            lattice_center[1].animate.set_value(0)
        )
        self.wait()
        self.play(av_lattice.animate.set_value(1))
        self.play(ah_lattice.animate.set_value(1))
        self.wait()
        self.move_camera(phi = 0*DEGREES)
        self.wait()
        self.play(
            ah_lattice.animate.set_value(1.3)
        )
        self.wait()
        self.play(
            av_lattice.animate.set_value(0.5)
        )
        self.wait()
        self.play(
            ah_lattice.animate.set_value(1),
            av_lattice.animate.set_value(1)
        )
        self.wait()
        self.play(
            lattice_center[1].animate.set_value(-3),
            recip_center[1].animate.set_value(3),
            recip_center[2].animate.set_value(0),
            lattice_center[2].animate.set_value(0),
        )
        self.wait()
        recip_copy = recip_2D.copy()
        for i in range(point_count_v):
            for j in range(point_count_h):
                recip_copy[i][j].set_opacity(0.5**(abs(i-int(point_count_v/2))+abs(j-int(point_count_h/2))))
        self.play(
            FadeTransform(recip_2D, recip_copy)
        )
        self.wait()
        self.play(
            FadeTransform(recip_copy, recip_2D)
        )
        self.wait()
        self.wait()

        lattice_1D.move_to(lattice_2D.get_center())
        recip_1D.move_to(recip_2D.get_center())

        lattice_1D_bracket = always_redraw(
            lambda: BraceBetweenPoints(lattice_1D[0].get_center(), lattice_1D[1].get_center()).set_color(BLUE)
        )
        recip_1D_bracket = always_redraw(
            lambda: BraceBetweenPoints(recip_1D[0].get_center(), recip_1D[1].get_center()).set_color(RED)
        )
        oneD_brackets = VGroup(lattice_1D_bracket, recip_1D_bracket)

        """
        lattice_constant = always_redraw(
            lambda: DecimalNumber(ah_lattice.get_value()).set_color(BLUE).move_to(lattice_1D_bracket.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        recip_constant = always_redraw(
            lambda: DecimalNumber(ah_lattice.get_value()).set_color(RED).move_to(recip_1D_bracket.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        """
        lattice_constant = always_redraw(
            lambda: MathTex("L").set_color(BLUE).move_to(lattice_1D_bracket.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        recip_constant = always_redraw(
            lambda: MathTex("2\\pi/L").set_color(RED).move_to(recip_1D_bracket.get_tip()+DOWN*0.2).rotate(angle=180*DEGREES, axis=OUT)
        )
        oneD_counters = VGroup(lattice_constant, recip_constant)


        self.play(
            ReplacementTransform(recip_2D, recip_1D),
            ReplacementTransform(lattice_2D, lattice_1D),
            ReplacementTransform(brackets, oneD_brackets),
            ReplacementTransform(twoD_counters, oneD_counters)
        )

        self.wait()
        self.add(real_1D_label := real_label.copy())
        self.add(recip_1D_label := recip_label.copy())
        self.remove(real_label)
        self.remove(recip_label)

        self.play(
            lattice_1D.animate.shift(DOWN*4).scale(0.7, about_point = [0,0,0]),
            real_1D_label.animate.shift(DOWN*4).scale(0.7, about_point = [0,0,0]),
            recip_1D.animate.shift(UP*4).scale(0.7, about_point = [0,0,0]),
            recip_1D_label.animate.shift(UP*4).scale(0.7, about_point = [0,0,0])
        )

        self.wait()

        statement = VGroup(
            Tex("Reciprocal space"),
            Tex("is the"),
            Tex("$\\textbf{Fourier Transform}$"),
            Tex("of real space")
        ).arrange(DOWN).scale(1.5).rotate(-90*DEGREES)
        statement.add(SurroundingRectangle(statement, color=GOLD, stroke_width=12, buff = 0.2))
        self.play(Write(statement))
        self.wait()
        fade_all(self)

class LatticeConstantDerivation(Scene):
    def construct(self):
        condition = Tex(
            "$\\mathbf{a}_i$",      #0
            "$\\mathbf{b}_j$",      #1
            " $=$ ",                   #2
            "$c\\delta$",            #3
            "$_i$",                   #4
            "$_j$"                   #5
        )
        VGroup(condition[0], condition[4]).set_color(BLUE)
        VGroup(condition[1], condition[5]).set_color(RED)

        lattice = get_2D_lattice(3, 3, i_ah = 1, i_av = 1)
        lattice_real = lattice['real'][0]
        lattice_recip = lattice['reciprocal'][0].set_opacity(0.3)

        self.play(DrawBorderThenFill(lattice_real))
        self.play(DrawBorderThenFill(lattice_recip))
        self.play(Write(condition.next_to(lattice_recip, UP, buff = 0.5)))

        ah = lattice['constants'][0]
        self.wait()
        self.play(ah.animate.set_value(1.2))
        self.wait()

class FourierDecomposition(Scene):
    def construct(self):
        srate = 500. #sampling rate in Hz

        #create arrays of frequencies, amplitudes, and phases to plot
        frex = np.array([3, 5, 15, 35])
        amplit = np.array([5, 9, 5, 7])
        phases = np.pi*np.array([1/7., 1/4., 1/2., -1/4.])
        time = [-1,1]
        decompositions = VGroup()

        sine_waves = np.zeros([len(frex), len(time)])
        for fi in range(len(frex)):
            # sine_waves[fi,:] = amplit[fi] * np.sin(2*np.pi*frex[fi]*time+phases[fi])

            ax = Axes(
                x_range = [time[0], time[-1]],
                x_length = 20
            ).scale(0.2)

            wave = ax.plot(
                lambda x: amplit[fi] * np.sin(2*np.pi*frex[fi]*x + phases[fi])
            ).set_color(RED)

            decompositions.add(wave)

        decompositions.arrange(DOWN, buff = 0.3).shift(LEFT*3)

        ax_sum = Axes(
            x_range = [time[0], time[-1]],
            x_length = 20
        ).scale(0.2).shift(RIGHT*2.5)

        def sum_func(x):
            val = 0
            for fi in range(len(frex)):
                val += amplit[fi] * np.sin(2*np.pi*frex[fi]*x + phases[fi])
            return val

        sum_plot = ax_sum.plot(
            sum_func
        ).set_color(BLUE)

        self.play(Create(decompositions), run_time = 2)
        self.wait()
        decomp_copies = VGroup()
        for decomposition in decompositions:
            decomp_copy = decomposition.copy()
            self.play(decomp_copy.animate.move_to(sum_plot.get_center()))
            decomp_copies.add(decomp_copy)
        self.wait()
        self.play(ReplacementTransform(decomp_copies, sum_plot))
        self.wait()
        self.play(ReplacementTransform(sum_plot.copy(), decompositions), run_time = 2)
        self.wait()
        self.play(
            Uncreate(VGroup(sum_plot, decomp_copies, decompositions), reverse = True),
            run_time = 1.5
        )
        self.wait()


class Introduction(Scene):
    def construct(self):

        gibbs_quote_1 = VGroup(
            Text("A mathematician may ", slant=ITALIC),
            Text("say anything he pleases,", slant=ITALIC),
            Text("but a physicist must be ", slant=ITALIC),
            Text("at least partially sane.", slant=ITALIC),
            Text("-J. Willard Gibbs").scale(0.9)
        ).scale(1.2).arrange(DOWN)

        gibbs_quote_2 = VGroup(
            Text("One of the principal objects of theoretical"),
            Text("research in my department of knowledge,"),
            Text("is to find the point of view from which"),
            Text("the subject appears in its greatest simplicity.")
        ).arrange(DOWN)

        self.play(Write(gibbs_quote_1), run_time=5)
        self.wait(2)
        fade_all(self)
        self.wait()

        #IMAGES BEFORE THIS WITH MRI AND SOLID STATE TEXTBOOK

        wth_is_kspace = VGroup(
            Text("What the heck is", font = "Helvetica"),
            Text("k-space?!", font = "Helvetica")
        ).scale(1.5).arrange(DOWN)
        self.play(DrawBorderThenFill(wth_is_kspace))
        reciprocal_statment = Text("reciprocal space?!", font = "Arial").scale(1.5).move_to(wth_is_kspace[1].get_center())
        self.wait()
        self.play(ReplacementTransform(wth_is_kspace[1], reciprocal_statment))

        self.wait()
        fade_all(self)
        """
        From here up until the next set of triple quotes. This code was taken from a StackExchange answer
        about creating piecharts: https://stackoverflow.com/questions/65783858/how-to-create-a-correct-pie-chart-with-manim
        """
        Sector.set_default(inner_radius=0, outer_radius=3, stroke_width=3, fill_opacity=0.7, stroke_color=GREY_BROWN)

        init_values = [20, 13]
        colors = [GREEN, BLUE]
        init_total = sum(init_values)
        init_angles = [360 * value / init_total for value in init_values]
        init_sangles = [sum(init_angles[:i]) for i in range(len(init_angles))]

        final_values = [90, 0.7]
        final_total = sum(final_values)
        final_angles = [360 * value / final_total for value in final_values]
        final_sangles = [sum(final_angles[:i]) for i in range(len(final_angles))]

        sectors = VGroup()
        for value, init_angle, init_sangle, final_angle, final_sangle, color in zip(init_values, init_angles,
                                                                                    init_sangles, final_angles,
                                                                                    final_sangles, colors):
            sector = Sector(
                start_angle=init_sangle * DEGREES,
                angle=init_angle * DEGREES,
                color=color,
            )
            sector.init_angle = init_angle
            sector.init_sangle = init_sangle
            sector.final_angle = final_angle
            sector.final_sangle = final_sangle
            sectors.add(sector)

        for sector in sectors:
            sector.save_state()

        def update_sector(sector, alpha):
            sector.restore()
            angle = interpolate(sector.init_angle * DEGREES, sector.final_angle * DEGREES, alpha)
            start_angle = interpolate(sector.init_sangle * DEGREES, sector.final_sangle * DEGREES, alpha)
            sector.become(
                Sector(
                    start_angle=start_angle,
                    angle=angle,
                    color=sector.color,
                )
            )

        """
        """

        chart_envelope = Circle(radius=3, color=GREY_BROWN, stroke_width=5)
        statement = VGroup(
            Tex("\\% of people"),
            Tex("that care"),
            Tex("about this")
        ).arrange(DOWN).next_to(chart_envelope, RIGHT, buff = 0.2).set_color(BLUE)

        # -------------------- Animation --------------------
        self.play(Write(statement))
        self.play(Create(sectors))
        self.add(chart_envelope)
        self.wait(1)
        self.play(
            UpdateFromAlphaFunc(sectors[0], update_sector),
            UpdateFromAlphaFunc(sectors[1], update_sector),
        )
        self.wait(1)
        fade_all(self)

        dont_understand_statement = VGroup(
            Text("I just dont", font = "Helvetica"),
            Text("get it :(", font = "Helvetica")
        ).scale(1.5).arrange(DOWN)
        self.play(Write(dont_understand_statement), run_time = 0.9)
        self.wait()
        do_understant_statement = Text("I get it! :)", font = "Helvetica").scale(1.5)
        self.play(ReplacementTransform(dont_understand_statement, do_understant_statement))
        self.wait()
        fade_all(self)

class ZoomScroll(ThreeDScene):
    def construct(self):
        lattice_2D = get_2D_lattice(9, 9, automatic = False)
        lattice_3D = VGroup()

        colors = [
            GRAY, GRAY_A, DARK_GRAY, GRAY_B, LIGHT_GRAY, GRAY, GRAY_A, GRAY_E, DARK_GRAY, LIGHT_GRAY, GRAY_C, GRAY_D, GRAY, GRAY_E,
            GRAY, GRAY_A, DARK_GRAY, LIGHTER_GRAY, LIGHTER_GREY, DARKER_GRAY, GRAY_BROWN, GRAY_E, GREY, GREY_B
        ]
        for i, color in enumerate(colors):
            lattice_3D.add(
                lattice_2D.copy().shift(IN*i).set_color(color)
            )

        self.play(Create(lattice_3D))
        self.wait()
        self.move_camera(
            frame_center = [0, 0, -8],
            zoom = 1.5,
            theta = 25*DEGREES,
            run_time = 12
        )
        self.wait()

class Drop(ThreeDScene):
    def construct(self):

        real_lattice = get_2D_lattice(4, 200, side_length = 1, i_ah = 5, i_av = 5, automatic = False).move_to([0,40,4])
        recip_lattice = get_2D_lattice(4, 200, side_length = 1, i_ah = 5, i_av = 5, automatic = False).set_color(RED).move_to([0,40,-4])

        self.add(real_lattice)
        self.add(recip_lattice)

        self.move_camera(
            phi = -90*DEGREES
        )

        self.move_camera(
            frame_center = [0, 100, 0],
            run_time = 2,
        )

        self.wait()

class BlochtoKspace(MovingCameraScene):
    def construct(self):
        bt_1D = VGroup(
            MathTex("\\textbf{Bloch's Theorem in 1D}"),     
            MathTex(
                "\\psi_k",        #0 
                "(",            #1
                "x",            #2
                "+",            #3
                "m",            #4
                "a",            #5
                ")",            #6
                "=",            #7
                "e^",           #8
                "{i",           #9
                "k",            #10
                "a}",           #11
                "\\psi_k",        #12
                "(",            #13
                "x",            #14
                ")"             #15
            ),
            MathTex(
                "\\psi_k",        #0
                "(",             #1
                "x",            #2
                ")",           #3
                "=",            #4
                "e^",           #5
                "{i",           #6
                "k",            #7
                "x}",           #8
                "u_k",            #9
                "(",            #10
                "x",            #11
                ")"             #12
            ),          
            MathTex(                                        
                "m = 0, \\pm 1, ... \\pm N/2",
                "\\text{ }",
                "k = \\frac{2\\pi m}{Na}"
            )
        ).arrange(DOWN)
        bt_1D.add(SurroundingRectangle(bt_1D, color = GOLD, buff = 0.3, stroke_width = 10))
        self.play(Write(bt_1D))
        self.wait()

        self.play(FadeIn(temp_rect1 := SurroundingRectangle(bt_1D[1], color=BLUE))),
        self.wait()
        self.play(FadeOut(temp_rect1))
        self.play(FadeIn(temp_rect2 := SurroundingRectangle(bt_1D[2], color=GREEN_B)))
        self.wait()
        self.play(FadeOut(temp_rect2))
        self.wait()
        self.play(FadeIn(temp_rect3 := SurroundingRectangle(bt_1D[3][0], color=LIGHT_BROWN)))
        self.wait()
        self.play(FadeOut(temp_rect3))
        self.wait()

        #NEW STUFF
        self.play(bt_1D[-2][-1].animate.set_color(GREEN))
        self.wait()
        not_equal = MathTex("\\neq \\frac{2\\pi m}{L} = k_m").next_to(bt_1D[-2][-1], RIGHT, buff = 0.3)
        self.play(
            FadeIn(
                not_equal_group := VGroup(
                    BackgroundRectangle(not_equal, buff = 0.2),
                    not_equal
                )
            )
        )
        self.wait()
        self.play(FadeOut(not_equal_group))
        self.wait()

        self.play(bt_1D.animate.shift(UP*2+LEFT*3.5).scale(0.7))
        self.wait()

        k_val_explanation = BulletedList(
            "Let $L = Na$: lattice length",
            "Then $k = \\frac{2\\pi m}{L}$",
            "Number of $k$ = N",
            "$k$ values spaced out by $\\frac{2\\pi}{L}$"
        ).next_to(bt_1D, DOWN, buff = 0.5).scale(0.9)

        for phrase in k_val_explanation:
            self.play(Write(phrase))
            self.wait()

        k_lattice = get_1D_lattice(11, i_ah = 10/11).rotate(135*DEGREES).move_to([2.5,0,0]).set_color(RED)
        k_lattice_braces = VGroup()
        k_lattice_labels = VGroup()
        for i in range(len(k_lattice)-1):
            k_lattice_braces.add(
                BraceBetweenPoints(
                    k_lattice[i].get_center(),
                    k_lattice[i+1].get_center()
                ).set_color(RED),
            )
            if i == int(len(k_lattice)/2):
                k_lattice_labels.add(
                    MathTex("\\delta k = 2\\pi/L").next_to(k_lattice_braces[-1], RIGHT, buff = 0).shift(UP*0.2).set_color(RED)
                )
            else:
                k_lattice_labels.add(
                    MathTex("\\delta k ").next_to(k_lattice_braces[-1], RIGHT, buff = 0).shift(UP*0.2).set_color(RED)
                )

        k_lattice_group = VGroup(k_lattice, k_lattice_braces, k_lattice_labels)
        self.play(DrawBorderThenFill(k_lattice_group))
        self.wait()
        self.play(
            *[
                FadeOut(mob) for mob in VGroup(k_lattice_group, k_val_explanation)
            ],
        )
        self.wait()

        periodicity_1D = VGroup(
            MathTex("\\psi_{\\left[k_m+\\frac{2\\pi m}{a}\\right]}(x + ma)"),
            MathTex("= e^{i(k_m+\\frac{2\\pi m}{a})a} \\psi_{\\left[k_m+\\frac{2\\pi m}{a}\\right]}(x)"),
            MathTex("=", "e^{i2\\pi m}", "e^{ik_ma} \\psi_{\\left[k_m+\\frac{2\\pi m}{a}\\right]}(x) "),
            Tex("$k_m$ can take on multiple"),
            Tex("values spaced $\\Delta k = \\frac{2\\pi}{a}$ apart")
        ).arrange(DOWN).next_to(bt_1D, DOWN, buff = -0.1).scale(0.9)
        self.play(Write(periodicity_1D[0]))
        self.wait()
        self.play(Write(periodicity_1D[1]))
        self.wait()
        self.play(Write(periodicity_1D[2]))
        self.wait()
        self.play(
            periodicity_1D[2][2].animate.shift(
                np.linalg.norm(periodicity_1D[2][1].get_left() - periodicity_1D[2][1].get_right())*LEFT
            ),
            FadeOut(periodicity_1D[2][1])
        )
        self.wait()
        self.play(Write(periodicity_1D[3:]))
        self.wait()

        periodicity_1D.add(
            VGroup(
                MathTex("\\textbf{Choose}"),
                MathTex("\\psi_{\\left[k_m+\\frac{2\\pi m}{a}\\right]}(x) = \\psi_{k_{m}}(x)}"),
                MathTex("u_{\\left[k_m+\\frac{2\\pi m}{a}\\right]}(x) = u_{k_{m}}(x)}")
            ).arrange(DOWN).move_to(bt_1D.get_center()+np.array([0,-0.1,0])).set_color(GOLD)
        )
        self.play(Write(periodicity_1D[-1]), FadeOut(bt_1D))
        self.wait()

        km_lattice = get_1D_lattice(5, i_ah = 3).rotate(135*DEGREES).move_to([2.5,0,0]).set_color(GREEN)
        km_lattice_braces = VGroup()
        km_lattice_labels = VGroup()
        for i in range(len(km_lattice)-1):
            km_lattice_braces.add(
                always_redraw(lambda i=i: BraceBetweenPoints(
                    km_lattice[i+1].get_center(),
                    km_lattice[i].get_center()
                ).set_color(GREEN)),
            )
            if i == int(len(k_lattice)/2):
                km_lattice_labels.add(
                    always_redraw(
                        lambda i=i: MathTex("\\Delta k_m = 2\\pi/a").next_to(km_lattice_braces[i].get_center(), LEFT, buff = 0.3).shift(DOWN*0.2).set_color(GREEN)
                    )
                )
            else:
                km_lattice_labels.add(
                    always_redraw(
                        lambda i=i: MathTex("\\Delta k_m").next_to(km_lattice_braces[i].get_center(), LEFT, buff = 0.3).shift(DOWN*0.2).set_color(GREEN)
                    )
                )

        km_lattice_group = VGroup(km_lattice, km_lattice_braces, km_lattice_labels)
        self.play(DrawBorderThenFill(km_lattice_group))
        self.wait()

        self.play(km_lattice.animate.scale(3))
        self.play(
            FadeIn(k_lattice_group),
        )
        self.wait()

        
        brillouin_zone = VGroup(
            MathTex("\\textbf{Brillouin}"),
            MathTex("\\textbf{Zone}")
        ).arrange(DOWN).move_to([4, 2.5, 0]).scale(0.9)
        gamma_circle = Circle(radius = 0.5).set_color(PINK).move_to(k_lattice[int(len(k_lattice)/2)].get_center())
        gamma_point = MathTex("k = 0").next_to(gamma_circle, DL, buff = 0).set_color(PINK)
        xp_circle = Circle(radius = 0.5).set_color(PURPLE_C).move_to(k_lattice[-1].get_center())
        xp_point = MathTex("k = \\frac{\\pi}{a}").next_to(xp_circle, RIGHT, buff = 0.3).set_color(PURPLE_C)
        xm_circle = Circle(radius = 0.5).set_color(PURPLE_C).move_to(k_lattice[0].get_center())
        xm_point = MathTex("k = \\frac{-\\pi}{a}").next_to(xm_circle, UP, buff = 0.2).set_color(PURPLE_C)

        recip_lattice_points = VGroup(
            gamma_circle, gamma_point, 
            BackgroundRectangle(xp_point, buff = 0.2), xp_circle, xp_point, 
            BackgroundRectangle(xm_point, buff = 0.2), xm_circle, xm_point 
        )
        self.play(Write(recip_lattice_points))
        self.wait()

        self.play(Write(brillouin_zone))
        self.wait()
        brillouin_zone.add(SurroundingRectangle(brillouin_zone))
        self.play(FadeIn(brillouin_zone[-1], buff = 0.2))
        self.wait()

        extended_lattice = k_lattice.copy().shift(LEFT*7.07 + UP*7.07).set_color(LOGO_RED)
        self.play(FadeIn(extended_lattice))
        self.wait()
        self.play(self.camera.frame.animate.set(width=30))
        self.wait()
        self.play(extended_lattice.animate(path_arc=PI/3).shift(DOWN*7.07+RIGHT*7.07), run_time = 2)
        self.wait()

        self.play(self.camera.frame.animate.set(width=14))
        self.wait()

        self.play(
            *[
                FadeOut(mob) for mob in VGroup(
                    extended_lattice, recip_lattice_points, km_lattice_group, k_lattice_group, periodicity_1D, brillouin_zone
                )
            ],
            bt_1D.animate.scale(1/0.7).move_to([0,0,0])
        )
        self.wait()
        """

        bt_3D = VGroup(
            MathTex("\\textbf{Bloch's Theorem in 3D}"),     
            Tex(
                "$\\psi_\\mathbf{k}$",        #0 
                "$($",            #1
                "$\\mathbf{r}$",  #2
                "$+$",            #3
                "$\\mathbf{R}$",  #4
                "$)$",            #5
                "$=$",            #6
                "$e$",            #7
                "$^i$",           #8
                "$^{\\mathbf{k}}$",#9
                "$^\\cdot$",         #10
                "$^{\\mathbf{a}}$",   #11
                "$\\psi_\\mathbf{k}$",        #12
                "$($",            #13
                "$\\mathbf{r}$",      #14
                "$)$"             #15
            ),
            Tex(
                "$\\psi_\\mathbf{k}$",        #0
                "$($",             #1
                "$\\mathbf{r}$",  #2
                "$)$",             #3
                "$=$",            #4
                "$e$",           #5
                "$^i$",           #6
                "$^\\mathbf{k}$", #7
                "$^\\cdot$",        #8
                "$^\\mathbf{r}$",   #9
                "$u_\\mathbf{k}$",            #10
                "$($",            #11
                "$\\mathbf{r}$",  #12
                "$)$"             #13
            ),          
            Tex(                                        
                "$\\mathbf{R}$ is any lattice vector",
                "\\\ $\\mathbf{k} = k_x\\hat{i} + k_y\\hat{j} + k_z\\hat{k}$"
            )
        ).arrange(DOWN)
        bt_3D.add(SurroundingRectangle(bt_3D, color = GOLD, buff = 0.3, stroke_width = 10))

        self.play(
            FadeOut(bt_1D),
            FadeIn(bt_3D)
        )
        self.wait()

        self.play(
            bt_3D.animate.move_to([-3.5, 0, 0]),
            bt_1D.animate.move_to([3.5, 0, 0]),
        )
        self.wait()

        equation_rects = VGroup()
        equation_rects.add(highlight(self, bt_1D[1][2], RED))
        equation_rects.add(highlight(self, bt_3D[1][2], RED))
        self.wait()
        equation_rects.add(highlight(self, bt_1D[1][4:6], BLUE))
        equation_rects.add(highlight(self, bt_3D[1][4], BLUE))
        self.wait()

        other_differences = VGroup(
            bt_1D[1][10:12], bt_3D[1][9:12], 
            bt_1D[1][-2], bt_3D[1][-2],
            bt_1D[2][2], bt_3D[2][2],
            bt_1D[2][7:9], bt_3D[2][7:10],
            bt_1D[2][11], bt_3D[2][12],
        )
        other_difference_colors = [ORANGE, PURPLE, YELLOW, GREEN, LOGO_RED]

        for i, color in enumerate(other_difference_colors):
            equation_rects.add(highlight(self, other_differences[2*i], color, rt = 0.2))
            equation_rects.add(highlight(self, other_differences[2*i+1], color, rt = 0.2))
        self.wait()
        equation_rects.add(
                highlight(self, bt_3D[3], PINK),
        )
        equation_rects.add(Rectangle(
            width = 6.8,
            height = 1.2,
            color = PINK
        ).move_to(bt_1D[3].get_center()))
        self.play(FadeIn(equation_rects[-1]))
        self.wait()

        """

class BlochReview(Scene):
    def construct(self):
        lattice = get_1D_lattice(13, i_ah = 1.5).shift(DOWN*2)
        self.add(lattice)

        wavefunction = MathTex(
            "\\psi(x + ma)",
            "=",
            "e^{ika}",
            "\\psi(x)"
        ).shift(UP*1 + LEFT*3.2).scale(1.2)
        wavefunction[2].set_color(RED)

        ax = NumberPlane(
            x_range = [-1.5, 1.5],
            y_range = [-1.5, 1.5],
            x_length = 4,
            y_length = 4,
        ).move_to([2.5, 1, 0])
        x_label = Tex("Re").next_to(ax.x_axis, RIGHT)
        y_label = Tex("Im").next_to(ax.y_axis, UP)

        vec = Arrow(
            start = [2.25, 1, 0],
            end = ax.c2p(1.2,0,0)
        ).set_color(RED)

        tip = always_redraw(
            lambda: Dot().scale(0.05).set_color(RED).move_to(vec.copy().scale(1.25).get_tip())
        )
        vec_trail = TracedPath(tip.get_center, dissipating_time=0.5, stroke_opacity = 0.3, stroke_width=4.2).set_color(RED)

        dot = Dot().scale(3).set_color(GREEN).shift(DOWN*2)
        trail = TracedPath(dot.get_center, dissipating_time=0.45, stroke_opacity = 0.38, stroke_width=15).set_color(GREEN)
        self.play(DrawBorderThenFill(dot))
        self.add(trail)
        self.wait()
        distances = np.array([1, 2, 1, -3, -1, -2, 2]) * 1.5

        self.play(FadeIn(wavefunction))
        self.wait()
        self.play(Create(VGroup(ax, x_label, y_label, vec)))
        self.wait()
        self.add(tip, vec_trail)

        for distance in distances:
            self.play(
                dot.animate.shift(RIGHT*distance),
                Rotate(vec, angle=distance, about_point=[2.5, 1, 0], axis=OUT),
                run_time = 1.5
            )

        self.wait()
        fade_all(self)
        self.wait()

        #K VALUE SHIFT VISUALIZATION
        ax = NumberPlane(
            x_range = [-1.5, 1.5],
            y_range = [-1.5, 1.5],
            x_length = 6,
            y_length = 6,
        )
        x_label = Tex("Re").next_to(ax.x_axis, RIGHT)
        y_label = Tex("Im").next_to(ax.y_axis, UP)

        first_angles = VGroup(
            MathTex("0").move_to([3.2, 0, 0]),
            MathTex("\\frac{\\pi}{2}").move_to([0, 3.2, 0]),
            MathTex("\\pi").move_to([-3.2, 0, 0]),
            MathTex("\\frac{3\\pi}{2}").move_to([0, -3.2, 0])
        ).scale(0.8)

        second_angles = VGroup(
            MathTex("0").move_to([3.2, 0, 0]),
            MathTex("\\frac{\\pi}{2}").move_to([0, 3.2, 0]),
            MathTex("\\pm\\pi").move_to([-3.2, 0, 0]),
            MathTex("-\\frac{\\pi}{2}").move_to([0, -3.2, 0])
        ).scale(0.8)

        tick_mark_backgrounds = VGroup()

        for i in range(len(first_angles)):
            tick_mark_backgrounds.add(
                BackgroundRectangle(first_angles[i], buff = 0.2)
            )


        vec = Arrow(
            start = [-0.25,0,0],
            end = ax.c2p(1.12, 0, 0)
        ).set_color(RED)

        tip = always_redraw(
            lambda: Dot().scale(0.05).set_color(RED).move_to(vec.get_end())
        )
        vec_trail = TracedPath(tip.get_center, dissipating_time=15, stroke_opacity = 0.85, stroke_width=5.5).set_color(RED)
        vec_2 = TracedPath(tip.get_center, dissipating_time=15, stroke_opacity = 0.85, stroke_width=5.5).set_color(RED)
        self.play(
            Write(VGroup(ax, x_label, y_label, vec, tip, tick_mark_backgrounds, first_angles))
        )
        self.add(vec_trail)

        self.play(Rotate(vec, angle = 2*PI, about_point = [0,0,0], axis = OUT), run_time = 6)
        self.wait()
        self.remove(vec_trail)
        self.wait()
        self.play(Rotate(vec, angle = PI, about_point = [0,0,0], axis = OUT), run_time = 2.5, rate_func = linear)
        self.wait()
        self.add(vec_2)
        self.play(ReplacementTransform(first_angles, second_angles))
        self.wait()
        self.play(Rotate(vec, angle = 2*PI, about_point = [0,0,0], axis = OUT), run_time = 6)
        self.wait()
        fade_all(self)




class ThreeDTorus(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi = 58*DEGREES,
        )

        donut = Torus(
            major_radius = 2.5,
            minor_radius = 1.25
        )
        self.add(donut)
        self.wait()

        big_circle = VGroup(
            Arc(
                start_angle = 0,
                angle = -PI,
                radius = 2.5 + 1.25,
                stroke_width = 10
            ).set_color(RED),
            Arc(
                start_angle = 0,
                angle = PI,
                radius = 2.5 + 1.25,
                stroke_width = 10,
                stroke_opacity = 0.3
            ).set_color(RED)
        )

        little_circle = Circle(
            radius = 1.25,
            stroke_width = 10,
            stroke_opacity = 0.9
        ).rotate(90*DEGREES, axis = UP).move_to([0,-2.5,0]).set_color(GREEN)


        self.play(Create(big_circle[0]), Create(big_circle[1]))
        self.wait()
        self.play(Create(little_circle))
        self.wait()

        self.play(
            donut.animate.rotate(angle = PI/6, axis = OUT),
            little_circle.animate.rotate_about_origin(angle = PI/6, axis = OUT),
            run_time = 3
        )
        self.wait()

class KSpaceMotivation(ThreeDScene):
    def construct(self):

        fourier_series = VGroup(
            Tex("If a function $f(x)$ is periodic with $L$"),
            Tex("then the Fourier Expansion of $f$ is: "),
            MathTex(
            "\\sum_{m=\\infty}^\\infty",    #0
            "A_m",                          #1
            "e^{ik_mx}, \\text{ }",         #2
            "k_m",                          #3      
            "=",                            #4
            "m",                            #5
            "\\frac{2\\pi}{L}"              #6
            ),
            Tex("$k_m$: Fourier component").set_color(BLUE),
            Tex("$A_m$: Fourier coefficient").set_color(PURPLE)
        ).arrange(DOWN)
        fourier_series[2][3].set_color(BLUE)
        fourier_series[2][1].set_color(PURPLE)

        self.play(Write(fourier_series[0]))
        self.wait()
        ax = Axes().shift(DOWN*1)
        func = ax.plot(
            lambda x: 2*np.sin(4*x) - np.cos(4*x) + 0.5*np.sin(8*x) + 0.3*np.cos(8*x) - 0.2*np.sin(12*x) + 0.4*np.cos(12*x) + 0.1*np.sin(16*x) - 0.3*np.cos(16*x)
        ).set_color(BLUE_D).shift(LEFT*0.2)
        self.play(Create(func))
        self.wait()
        dot = Dot().move_to(ax.c2p(-0.8, 2*np.sin(4*(-0.8)), 0))
        distances = [1, 2, 1, 0, -1, -2]
        copies = VGroup(dot)
        for distance in distances:
            copies.add(dot.copy().set_opacity(0.4))
            self.add(copies[-1])
            points = [ax.c2p(-0.8 + distance*np.pi/2-PI/2, 2*np.sin(4*(-0.8 + distance*PI)), 0), ax.c2p(-0.8 + distance*np.pi/2, 2*np.sin(4*(-0.8 + distance*PI)), 0)]
            distance_brace = BraceBetweenPoints(dot.copy().get_center(), dot.copy().move_to(ax.c2p(-0.8 + distance*np.pi/2, 2*np.sin(4*(-0.8 + distance*PI)), 0)).get_center())
            L = MathTex("L").move_to(distance_brace.get_center() + np.array([0,1,0])*copysign(0.4, points[1][0]-points[0][0]))
            self.play(
                FadeIn(distance_brace),
                FadeIn(L),
                dot.animate.move_to(ax.c2p(-0.8 + distance*np.pi/2, 2*np.sin(4*(-0.8 + distance*PI)), 0)), 
            )
            self.play(FadeOut(distance_brace), FadeOut(L))

        copies.add(dot.copy().set_opacity(0.4))
        self.add(copies[-1])
        distance_braces = VGroup()
        for distance in range(-1, 4):
            points = [ax.c2p(-0.8 + distance*np.pi/2-PI/2, 2*np.sin(4*(-0.8 + distance*PI)), 0), ax.c2p(-0.8 + distance*np.pi/2, 2*np.sin(4*(-0.8 + distance*PI)), 0)]
            brace = BraceBetweenPoints(points[0], points[1])
            distance_braces.add(
                brace,
                MathTex("L").move_to(brace.get_center() + np.array([0, -copysign(0.4, points[1][0]-points[0][0]), 0]))
            )

        self.play(
            Create(distance_braces),
            dot.animate.move_to(ax.c2p(-0.8 + 3*np.pi/2, 2*np.sin(4*(-0.8 + 3*PI)), 0)), 
            run_time = 2
        )
        self.wait()

        self.play(Unwrite(VGroup(func, copies, distance_braces)))
        self.wait()

        for phrase in fourier_series[1:]:
            self.play(Write(phrase))
            self.wait()
        self.play(
            *[
                mob.animate.shift(UP*2+LEFT*4.5).scale(0.4, about_point=[-4.5,2,0]) for mob in VGroup(
                    fourier_series[0:2], fourier_series[2][0:3], fourier_series[3:]
                )
            ],
            fourier_series[2][3:].animate.move_to([0,0,0])
        )
        self.wait()
        fourier_series.add(k_basis_rect := SurroundingRectangle(fourier_series[2][6], color = RED))
        self.play(FadeIn(k_basis_rect))
        self.wait()

        k_basis = VGroup(
            k_basis_rect,
            fourier_series[2][3:],
        )   

        lattice = get_2D_lattice(5, 5)
        real_lattice = lattice['real'][0]
        recip_lattice = lattice['reciprocal'][0]

        middle_line = VGroup(VGroup(recip_lattice[0][-3], recip_lattice[1][-3], recip_lattice[2][-3], recip_lattice[3][-3], recip_lattice[4][-3])).copy()
        self.play(
            DrawBorderThenFill(middle_line),
            k_basis.animate.next_to(recip_lattice, UP, buff = 0.1)
        )
        self.wait()

        recip_lattice_brace = BraceBetweenPoints(recip_lattice[-3][-3].get_center(), recip_lattice[-4][-3].get_center()).set_color(RED)

        recip_lattice_counter = always_redraw(
            lambda: MathTex(
                "\\delta k",
                "=",
                "\\frac{2\\pi}{L}",
                
            ).next_to(recip_lattice_brace.get_center(), RIGHT).set_color(recip_lattice_brace.get_color())
        )

        self.play(
            Create(recip_lattice_counter), Create(recip_lattice_brace)
        )
        self.wait()
        self.play(
            Create(recip_lattice),
            recip_lattice_brace.animate.become(
                BraceBetweenPoints(recip_lattice[-3][-1].get_center(), recip_lattice[-4][-1].get_center()).set_color(RED)
            )
        )
        self.remove(middle_line)
        self.wait()

        reciprocal_statement = VGroup(
            Tex("Reciprocal Space/Lattice"),
            Tex("(Fourier Space)").scale(0.7)
        ).scale(1.3).arrange(DOWN)
        reciprocal_lattice = VGroup()
        reciprocal_lattice.add(BackgroundRectangle(reciprocal_statement, buff = 0.3))
        reciprocal_lattice.add(SurroundingRectangle(reciprocal_statement, color = GOLD, buff = 0.3))
        reciprocal_lattice.add(reciprocal_statement)
        self.play(Write(reciprocal_lattice))
        self.wait(2)
        
        self.play(
            *[
                Unwrite(mob) for mob in VGroup(
                    recip_lattice,
                    reciprocal_lattice,
                    recip_lattice_counter,
                    recip_lattice_brace,
                    fourier_series,

                )
            ],
            run_time = 1
        )

        self.wait()

        """
        lattice_dict = get_2D_lattice(200, 100, automatic=False, shape = Square)

        real_lattice = lattice_dict

        self.add(real_lattice)
        self.wait()
        self.move_camera(
            zoom = 0.1
        )
        self.wait()
        """

class PlaneWaves(ThreeDScene):
    def construct(self):
        axis = Axes(
            x_range = [-12, 12],
            y_range = [-6, 6],
            x_length = 12,
            y_length = 4
        )

        amp = ValueTracker(2)
        k = ValueTracker(1)
        phi = ValueTracker(0)

        wave = always_redraw(lambda: axis.plot(
            lambda x: amp.get_value() * np.cos(k.get_value()*x + phi.get_value())
        ).set_color(BLUE))

        self.play(
            Create(axis)
        )
        self.wait()
        self.play(Create(wave))
        self.wait()

        equation = MathTex(
            "f(x)",     #0
            "=",        #1
            "A",        #2    
            "\\cos(",   #3
            "k",        #4
            "x",        #5
            "+",        #6
            "\\phi",    #7
            ")"         #8
        ).shift(UP*2.5)
        k_equation = MathTex("k = 2\\pi /L").shift(UP*1.8)

        equation_group = VGroup(BackgroundRectangle(equation), equation)

        self.play(Write(equation_group))
        self.wait()

        #Begin showing parameters of plane wave
        equation_group.add(SurroundingRectangle(equation[2], PURPLE))
        self.play(FadeIn(equation_group[-1]))
        self.wait()
        self.play(
            amp.animate.set_value(3)
        )
        self.play(
            amp.animate.set_value(1)
        )
        self.play(
            amp.animate.set_value(2)
        )
        self.wait()
        self.play(FadeOut(equation_group[-1]))

        equation_group.add(SurroundingRectangle(equation[4], PURPLE))
        self.play(FadeIn(equation_group[-1]))
        self.wait()
        self.play(
            k.animate.set_value(2)
        )
        self.play(
            k.animate.set_value(3)
        )
        self.play(
            k.animate.set_value(1)
        )
        self.wait()
        self.play(FadeOut(equation_group[-1]))
        
        equation_group.add(SurroundingRectangle(equation[7], PURPLE))
        self.play(FadeIn(equation_group[-1]))
        self.wait()
        self.play(
            phi.animate.set_value(np.pi/2)
        )
        self.play(
            phi.animate.set_value(PI)
        )
        self.play(
            phi.animate.set_value(0)
        )
        self.wait()
        self.play(FadeOut(equation_group[-1]))

        self.wait()

        self.play(Write(k_rect := BackgroundRectangle(k_equation)), Write(k_equation))
        self.wait()
        
        #COMPLEX PART
        complex_equation = MathTex("f(x) = Ae^{ikx}").move_to(equation.get_center())
        axis_3D = ThreeDAxes(
            x_range = [-12, 12],
            y_range = [-6, 6],
            z_range = [-12, 12],
            x_length = 12,
            y_length = 4,
            z_length = 12
        )

        complex_wave = ParametricFunction(
            lambda t: axis_3D.c2p(
                t,
                2*np.cos(t),
                2*np.sin(t)
            ),
            t_range = [-4*np.pi, 4*np.pi],
        ).set_color(BLUE)

        complex_stuff = VGroup(
            axis_3D,
            complex_wave,
        )

        self.play(ReplacementTransform(equation, complex_equation))
        self.wait()
        self.play(DrawBorderThenFill(axis_3D))
        self.wait()
        self.play(ReplacementTransform(wave, complex_wave))
        self.remove(axis)
        self.remove(wave)
        self.play(
            *[
                mob.animate.rotate(
                    axis = UP,
                    angle = PI/5,
                    about_point = [0,0,0]
                ) for mob in complex_stuff
            ],
            run_time = 5
        )
        self.wait()
        self.play(
            *[
                mob.animate.rotate(
                    axis = UP,
                    angle = -PI/5,
                    about_point = [0,0,0]
                ) for mob in complex_stuff
            ],
            run_time = 5
        )
        self.wait()



class DimensionalGeneralization(Scene):
    def construct(self):
        fourier_1D = VGroup(
            Tex("$\\textbf{1D Fourier Series}$"),
            MathTex(
                "f(x) = "
                "\\sum_{m=\\infty}^\\infty",    #0
                "A_m",                          #1
                "e^{ik_mx}, \\text{ }",         #2
                "k_m",                          #3      
                "=",                            #4
                "2\\pi",                            #5
                "m /L"              #6
            ),
        ).arrange(DOWN).shift(UP*1.5)

        fourier_3D = VGroup(
            Tex("$\\textbf{3D Fourier Series}$"),
            Tex("hi").set_opacity(0).scale(0.5),
            Tex(
                "$$f(x) = \\sum_{m} A_m e^{i\\mathbf{k}_m \\cdot \\mathbf{r}}$$"
            ),
            Tex("$$\\mathbf{k}_m = k_{xm}\\hat{i} + k_{ym}\\hat{j} + k_{zm}\\hat{k}$$")
        ).arrange(DOWN).shift(DOWN*1.5)

        self.play(Write(fourier_1D))
        self.wait()
        self.play(Write(fourier_3D))
        self.wait()

        twoD_requirements = VGroup(
            Tex("If some function $f(\\mathbf{r})$ is $\\textbf{periodic}$ with $\\mathbf{R}$"),
            Tex("Where $\\mathbf{R} = n_1\\mathbf{a} + n_2\\mathbf{b} + n_3\\mathbf{c}, \\text{ } (n_1, n_2, n_3)\\in \\mathbb{Z}$"),
            Tex("then the 3D Fourier Series of $f$ is:"),
            VGroup(
                Tex("$$f(x) = \\sum_{m} A_m e^{i\\mathbf{k}_m \\cdot \\mathbf{r}} \\Leftrightarrow $$"),
                Tex("$$ A_m = \\frac{1}{\\nu}\\int_\\nu f(\\mathbf{r}) e^{-i \\mathbf{k}_m \\cdot \\mathbf{r}} d\\mathbf{r}$$").set_color(PURPLE_A),
            ).arrange(RIGHT),
            Tex("$$\\mathbf{k}_m = k_{xm}\\hat{i} + k_{ym}\\hat{j} + k_{zm}\\hat{k}, \\text{ } m \\in \\mathbb{Z}$$"),
            VGroup(
                Tex("For the reciprocal lattice: "),
                Tex("$$e^{i \\mathbf{k}_m \\cdot \\mathbf{R}} = 1$$")
            ).arrange(RIGHT)
            
        ).arrange(DOWN)

        self.play(FadeOut(fourier_1D))
        self.wait()
        self.play(
            Unwrite(fourier_3D, reverse = True),
            Write(twoD_requirements[0])
        )
        self.wait()
        for requirement in twoD_requirements[1:]:
            self.play(Write(requirement))
            self.wait()



class SchrodingerEquation(Scene):
    def construct(self):
        schrodinger_eqn = VGroup(
            Tex(
                "Free Electron ",
                "Schrodinger Equation"
            ),
            MathTex(
                "\\left[",                  #0
                "\\frac{-\\hbar^2}{2m",     #1
                "^*}",                      #2
                "\\frac{\\partial^2}{\\partial x^2}", #3
                "+",                        #4
                "V(x)",                     #5
                "\\right]",                 #6
                "\\psi_k(x)",              #7
                "=",                        #8
                "E",                        #9
                "\\psi_k(x)"                #10
            ),
            Tex("Calculate $E$ for all $k$"),
            Tex("Assume $V(x)$ is constant")
        ).arrange(DOWN, buff = 1)
        schrodinger_eqn[1][2].set_opacity(0)
        k_schrodinger_equation = Tex("$$\\hat{H}(k)u_k(x) = E(k)u_k(x)$$")
        for i, part in enumerate(schrodinger_eqn[:len(schrodinger_eqn)-1]):
            if i == 1:
                self.play(
                    Write(VGroup(part[0:2], part[3:]))
                )
                self.wait()
            else:
                self.play(Write(part))
                self.wait()
        self.wait()
        self.play(Write(schrodinger_eqn[3]))
        self.wait()

        self.wait()
        self.play(VGroup(schrodinger_eqn).animate.move_to([-3,1,0]).scale(0.7))
        
        self.wait()


        ax = Axes(
            x_range = [-3, 3],
            y_range = [-0.5, 4],
            x_length = 5,
            y_length = 6
        ).move_to([2.5,0,0])
        x_label = Tex("$k$").next_to(ax.x_axis, RIGHT)
        y_label = Tex("$E(k)$").next_to(ax.y_axis, UP)

        self.play(Create(VGroup(ax, x_label, y_label)))
        self.wait()

        free_electron_Ek = ax.plot(
            lambda x: 1.3*x**2 + 1,
            x_range = [-1.2, 1.2]
        ).set_color(BLUE)
        self.play(Create(free_electron_Ek))
        self.wait()
        free_electron_func = VGroup(
            MathTex("E(k) = V + \\frac{-\\hbar^2k^2}{2m}").move_to(ax.c2p(0, 3.4, 0)),
        )
        free_electron_func.add(BackgroundRectangle(free_electron_func))
        free_electron_func = VGroup(
            free_electron_func[1],
            free_electron_func[0].set_color(BLUE)
        ).scale(0.7)

        self.play(Write(free_electron_func))
        self.wait()

        #EFFECTIVE MASS
        crystal_electron = Tex("Crystal Electron").move_to(schrodinger_eqn[0][0].get_center() + DOWN*0.03 + LEFT*0.25).scale(0.7)
        free_to_crystal_arrow = Arrow(start = free_electron_func.get_left(), end = schrodinger_eqn[0][0].get_right()).set_color(RED)
        self.play(FadeIn(free_to_crystal_arrow))
        self.wait()
        self.play(FadeTransform(schrodinger_eqn[0][0], crystal_electron), FadeOut(free_to_crystal_arrow))
        self.wait()

        k_schrodinger_equation = k_schrodinger_equation.copy()
        self.play(FadeTransform(schrodinger_eqn[1], k_schrodinger_equation.scale(0.7).move_to(schrodinger_eqn[1].get_center())))
        schrodinger_eqn.remove(schrodinger_eqn[1])
        schrodinger_eqn.add(k_schrodinger_equation)
        self.wait()
        potential_fit = Tex("Fit $V(x)$ to data").move_to(schrodinger_eqn[2].get_center()).scale(0.7)
        self.play(FadeTransform(schrodinger_eqn[2], potential_fit))
        schrodinger_eqn.add(potential_fit)
        self.wait()
        # self.play(schrodinger_eqn[1][2].animate.set_opacity(1).set_color(PURE_GREEN))
        # self.play(Circumscribe(schrodinger_eqn[1][2]))
        self.wait()

        crystal_electron_Ek = ax.plot(
            lambda x: -1.5*np.cos(PI/2 * x)**2 + 2.5,
        ).set_color(GREEN)

        self.play(Create(crystal_electron_Ek))
        self.wait()
        self.play(free_electron_func.animate.shift(LEFT*1.5))

        crystal_electron_func = VGroup(
            Tex("$E(k)$ is periodic").move_to(ax.c2p(0, 3.4, 0)),
        )
        crystal_electron_func.add(BackgroundRectangle(crystal_electron_func))
        crystal_electron_func = VGroup(
            crystal_electron_func[1],
            crystal_electron_func[0].set_color(GREEN)
        ).scale(0.7).next_to(free_electron_func, RIGHT, buff=0.5)
        self.play(Write(crystal_electron_func))
        self.wait()

        labels = VGroup(
            MathTex("\\pi/a").move_to(ax.c2p(1, -0.35, 0)),
            MathTex("-\\pi/a").move_to(ax.c2p(-1, -0.35, 0))
        )
        pi_over_a_lines = VGroup(
            DashedLine(
                start = ax.c2p(1,0,0), end = ax.c2p(1, 3, 0)
            ),
            DashedLine(
                start = ax.c2p(-1,0,0), end = ax.c2p(-1, 3, 0)
            )
        ).set_opacity(0.6)

        self.play(Write(labels))
        self.play(Write(pi_over_a_lines))
        self.wait()

        BZ_brace = BraceBetweenPoints(
            ax.c2p(-1, 1.2, 0),
            ax.c2p(1, 1.2, 0)
        ).set_color(GOLD)
        BZ_statement = Tex("Brillouin Zone").move_to([-3, -2, 0]).set_color(GOLD)
        pointer = Line(
            start = BZ_brace.get_tip(),
            end = BZ_statement.get_right() + RIGHT*0.1
        ).set_color(GOLD)

        self.play(Create(BZ_brace))
        self.wait()
        self.play(DrawBorderThenFill(pointer))
        self.wait()
        self.play(Write(BZ_statement))
        self.wait()

        approx_parabola = ax.plot(
            lambda x: 1 + 3.7011*x**2,
            x_range = [-0.7,0.7]
        ).set_color(LIGHT_BROWN)
        self.play(Create(approx_parabola))
        self.wait()

        approx_energy = MathTex("E(x) \\approx V + \\frac{-\\hbar^2 k^2}{2","m_{\\text{eff}}}").set_color(GOLD).move_to(crystal_electron_func.get_center()).scale(0.7)
        self.play(FadeTransform(crystal_electron_func, approx_energy))
        self.wait()
        self.play(approx_energy[1].animate.set_color(PURE_RED))
        self.wait()

class RealvsReciprocal(Scene):
    def construct(self):
        takeaway_1 = VGroup(
            Tex("Real Space: Period", " $N$"),
            Tex("k-Space: Period"," $2\\pi /a$")
        ).scale(2).arrange(DOWN)
        takeaway_1[0][1].set_color(BLUE)
        takeaway_1[1][1].set_color(RED)
        self.play(Write(takeaway_1))
        self.wait()

        comparison = MathTex("N", ">>", "\\frac{2\\pi}{a}").scale(3)
        comparison[0].set_color(BLUE)
        comparison[2].set_color(RED)

        self.play(ReplacementTransform(takeaway_1, comparison))
        self.wait()

        self.wait()

class Summary(ThreeDScene):
    def construct(self):
        summary = BulletedList(
            "Fourier components form a '$k$-space' lattice",
            "Reciprocal space is the FT of real space",
            "$\\psi_k$ periodic in real space and $k$ space (by choice)",
            "Brillouin Zone exists from $[\\frac{-\\pi}{a},\\frac{\\pi}{a}]$ in $k$-space",
            "All $\\psi_k$ are fully described in Brillouin Zone",
            "Periodic wave-function $\\implies E(k)$ periodic in $k$-space"
        )
        self.camera.frame_center = np.array([0, summary[0].get_center()[1]-1, 0])

        for item in summary:
            self.move_camera(
                frame_center = np.array([0, item.get_center()[1]-1.2, 0]),
            )
            self.play(Write(item))
            self.play(item.animate.scale(1.15).set_color(GOLD))
            self.wait()
            self.play(item.animate.scale(1/1.15).set_color(WHITE))

        self.wait()
        self.move_camera(
            frame_center = np.array([0, item.get_center()[1]+1, 0]),
        )
        self.wait()

class DeltaFunctions(Scene):
    def construct(self):

        lattice = get_1D_lattice(11, i_ah = 1.2)

        ax1 = Axes(
            x_range = [-5, 5],
            y_range = [-3, 3],
            x_length = 12,
            y_length = 3
        )

        self.play(DrawBorderThenFill(lattice))
        self.wait()

        self.play(DrawBorderThenFill(ax1))
        self.wait()

        deltas = VGroup()
        for point in lattice:
            deltas.add(Line(
                start = point.get_center(),
                end = point.get_center() + np.array([0,7,0])
            ).set_color(GREEN).set_opacity(0.5))

        self.play(LaggedStart(*[
            Create(delta) for delta in deltas
        ]))
        self.wait()
        cutoff_rect = Rectangle(
            width = 20
        ).set_fill(BLACK).set_color(BLACK).set_opacity(1).shift(UP*7)
        self.play(cutoff_rect.animate.move_to([0, 3, 0]))
        self.wait()

        real_space_group = VGroup(
            ax1, lattice, deltas, cutoff_rect
        )

        func_representation = MathTex("f(x) = \\sum_{m=-\\infty}^{\\infty} \\delta(x - ma) ").set_color(BLUE)
        self.play(
            real_space_group.animate.shift(UP*3).scale(0.7, about_point = [0,0,0]),
            Write(func_representation)
        )
        self.wait()

        fourier_transform_generic = MathTex("F\\{f(x)\\} = \\int_{-\\infty}^{\\infty}f(x)e^{-ikx}dx").next_to(func_representation, DOWN)
        FT_specific = VGroup(
             MathTex("F\\{f(x)\\} = \\int_{-\\infty}^{\\infty}\\left[\\sum_{m=-\\infty}^{\\infty} \\delta(x - ma)\\right]e^{-ikx}dx").next_to(func_representation, DOWN),
             MathTex("F\\{f(x)\\} = \\sum_{m=-\\infty}^{\\infty} \\left[ \\int_{-\\infty}^{\\infty}\\delta(x - ma)e^{-ikx}dx\\right]").next_to(func_representation, DOWN),
             MathTex("F\\{f(x)\\} = \\sum_{m=-\\infty}^{\\infty} e^{-ikma}").next_to(func_representation, DOWN),
        )
        dirac_comb = MathTex("F\\{f(x)\\} = \\frac{2\\pi}{a}\\sum_{m=-\\infty}^{\\infty} \\delta\\left(k -", "\\frac{2\\pi m}{a}","\\right)").next_to(func_representation, DOWN)
        self.play(Write(fourier_transform_generic))
        self.wait()
        self.play(FadeTransformPieces(fourier_transform_generic, FT_specific[0]))
        self.wait()
        self.play(FadeTransformPieces(FT_specific[0], FT_specific[1]))
        self.wait()
        self.play(FadeTransformPieces(FT_specific[1], FT_specific[2]))
        self.wait()
        self.play(FadeTransform(FT_specific[2], dirac_comb))
        self.wait()
        self.play(Create(temp_rect := SurroundingRectangle(dirac_comb[1], buff = 0.1, color = GOLD)))
        self.wait()
        self.play(
            VGroup(func_representation, dirac_comb, temp_rect).animate.scale(0.7, about_point = [0,0,0]).shift(UP*0.8),
        )
        self.wait()
        ax2 = ax1.copy()

        recip = lattice.copy().set_color(RED)

        recip_deltas = VGroup(
            Line(
                start = [-10,0,0],
                end = [10,0,0]
            ).set_color(GREEN).set_opacity(0.5)
        )
        for point in recip:
            recip_deltas.add(Line(
                start = point.get_center(),
                end = point.get_center() + np.array([0,1.4,0])
            ).set_color(GREEN).set_opacity(0.5))

        recip_space_group = VGroup(
            ax2, recip, recip_deltas
        ).move_to([0,-3,0])

        self.play(
            Write(recip_space_group)
        )
        self.wait()

class FinalProblems(MovingCameraScene):
    def construct(self):
        problems = VGroup(
            Tex("Problems to Ponder").set_color(GOLD).scale(1.5),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            #PROBLEM 1
            Tex("I said previously that shifting from $s = 0,...N-1$"),
            Tex("to $m = -N/2...N/2$ changes nothing. In fact, we have"),
            Tex("to omit one value from $m$ in order to keep the same range"),
            Tex("I informally left this out, but why doesn't it"),
            Tex("change very much from a practical sense?"),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("The wavefunctions $\\psi_k$ for each value of $k$"),
            Tex("all solve the Schrodinger equation. Show that"),
            Tex("a sum over all $\\psi_k$ forms a Fourier Series"),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("For a free electron, the value $\\hbar k = p$ is the momentum."),
            Tex("While this does not hold for the crystal electron"),
            Tex("Show that $\\hbar k$ as a quantity is still conserved"),
            Tex("under a scattering process that changes $\\hbar k$ by"),
            Tex("an integer multiple of $\\hbar 2\\pi/a$, this is referred"),
            Tex("to as translation by a reciprocal lattice vector"),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("hi").set_opacity(0),
            Tex("That's all for now. Bye bye! :)")
        ).arrange(DOWN)
        self.play(FadeIn(problems))
        self.camera.frame_center = problems[0].get_center()
        self.wait()
        self.play(
            self.camera.frame.animate.shift(DOWN*abs(np.array(problems.get_top() - problems.get_bottom())[1])),
            run_time = 90,
            rate_func = linear
        )
        self.wait()

        

