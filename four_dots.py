from manim import *
import numpy as np
import random as rd

#FADE OUT EVERYTHING ON SCREEN
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

equilateral_array = np.array([
        [-0.5,np.sqrt(3)/4,0],
        [0.5,np.sqrt(3)/4,0],
        [0,-np.sqrt(3)/4,0],
    ]
)

rhombus_array = np.array([
        [0,1,0],
        [-1/np.sqrt(3),0,0],
        [1/np.sqrt(3),0,0],
        [0,-1,0]
    ]
)

square_array = np.array([
        [-0.5,0.5,0],
        [-0.5,-0.5,0],
        [0.5,-0.5,0],
        [0.5,0.5,0]
    ]
)

centered_triangle_array = np.array([
        [0,0,0],
        [0,2*np.sqrt(3)/6,0],
        [-0.5,-np.sqrt(3)/6,0],
        [0.5,-np.sqrt(3)/6,0]
    ]
)

kite_array = np.array([
        [0,0.5,0],
        [0,-0.5,0],
        [-0.5,0.5-(1-np.sqrt(3)/2),0],
        [0.5,0.5-(1-np.sqrt(3)/2),0],
    ]
)

isosceles_array = 0.75*np.array([
        [0,0.5,0],
        [0,-0.5,0],
        -np.array([-0.5,-(0.5-(1+np.sqrt(3)/2)),0]),
        -np.array([0.5,-(0.5-(1+np.sqrt(3)/2)),0]),
    ]
)

trapezoid_array = 0.8*np.array([
        [-0.5, 0.5*np.sin(2*PI/5), 0],
        [0.5, 0.5*np.sin(2*PI/5), 0],
        [-0.5-np.cos(2*PI/5), -0.5*np.sin(2*PI/5), 0],
        [0.5+np.cos(2*PI/5), -0.5*np.sin(2*PI/5), 0]
    ]
)

#rotation matrix goes hard ngl
pentagon_array = np.empty(shape=(5,3))
for i in range(5):
    pentagon_array[i] = np.matmul([[np.cos(2*i*np.pi/5),-np.sin(2*i*np.pi/5),0],
                                   [np.sin(2*i*np.pi/5),np.cos(2*i*np.pi/5),0],
                                   [0,                0,                    0]],  
                                   [0,1,0])

class Homework_vs_Research(Scene):
    def construct(self):
        quote = Text("Geometry is not true; it is advantageous.")
        poincare_citation = Text("-Henri Poincaré", slant=ITALIC).scale(0.8).shift(DOWN*0.7)

        quote = VGroup(quote, poincare_citation)
        self.play(Write(quote), run_time = 4)
        self.wait(3)
        fade_all(self)
        self.wait(3)

        homework_v_research = Tex(
            "Homework", #0
            "Research", #1
            "Puzzles",  #2
        )

        self.play(Write(homework_v_research[0].move_to((0,0,0))))
        self.wait()
        
        self.play(homework_v_research[0].animate.move_to([-4,2.5,0]))
        self.wait()

        self.play(Write(homework_v_research[1].move_to([3.8,2.5,0])))
        self.wait()


        divider = Line(start=[0,4,0],end=[0,-4,0])
        self.play(Create(divider))
        self.wait(2)

        similarities = Text(
            "Critical Thinking ",    #0
            "Difficult ",            #1
        ).shift(UP*2)
        differences = Text(
            "Known Answer ", "Murky Answers ", #0, 1
            "1 Textbook ", "Scientific Papers, More Textbooks ", #2, 3
        ).move_to(similarities.get_center()+DOWN*1)

        self.play(Create(similarities[0]), Create(similarities[0].move_to([2,1,0])))
        self.wait()
        self.play(Create(differences))
        self.wait()

class ThreeTypesOfProblems(Scene):
    def construct(self):
        quote = Text("Geometry is not true; it is advantageous.")
        poincare_citation = Text("-Henri Poincaré", slant=ITALIC).scale(0.8).shift(DOWN*0.7)

        quote = VGroup(quote, poincare_citation)
        self.play(Write(quote), run_time = 4)
        self.wait(3)
        fade_all(self)
        self.wait()

        axis = Axes(
            x_range = [0,12],
            y_range = [0,8],
            x_length = 10,
            y_length = 6,
        )
        axis_center = axis.get_center()
        dividers = VGroup(
            DashedLine(start=axis_center+[0,axis.y_length/2,0], end=axis_center-[0,axis.y_length/2,0]),
            DashedLine(start=axis_center+RIGHT*axis.x_length/2, end=axis_center+LEFT*axis.x_length/2)
        ).set_color(PURPLE_E)
        
        self.play(Create(axis))
        axis_labels = axis.get_axis_labels(Text("Straightforwardness").scale(0.5), Text("Difficulty").scale(0.5))
        self.play(FadeIn(axis_labels[0].set_color(LOGO_WHITE)), FadeIn(axis_labels[1].set_color(LOGO_WHITE)))
        self.wait()

        #FERMAT'S LAST THEOREM
        fermat_point = Dot(axis.c2p(0.4,7.6)).set_color(RED)
        fermat_theorem = Text("Fermat's Last Theorem").scale(0.5).next_to(fermat_point, DOWN, buff=0.1).set_color(RED)
        self.play(FadeIn(fermat_point))
        self.play(Write(fermat_theorem))    
        fermat = VGroup(fermat_point, fermat_theorem) 

        #ARITHMETIC
        arithmetic_point = Dot(axis.c2p(11,0.2)).set_color(GREEN)
        arithmetic_theorem = Text("Arithmetic").scale(0.5).next_to(arithmetic_point, UL, buff=0.05).set_color(GREEN)
        self.play(FadeIn(arithmetic_point))
        self.play(Write(arithmetic_theorem))     
        arithmetic = VGroup(arithmetic_point, arithmetic_theorem)

        #TRAVELING SALESMAN 
        ts_point = Dot(axis.c2p(6,5))
        ts_theorem = Text("Traveling Salesman").scale(0.5).next_to(ts_point, UL, buff=0.05)
        ts = VGroup(ts_point, ts_theorem)

        #MONTY HALL
        monty_hall_point = Dot(axis.c2p(6,2))
        monty_hall_theorem = Text("Monty Hall Problem").scale(0.5).next_to(monty_hall_point, DL, buff=-0.05)
        monty_hall = VGroup(monty_hall_point, monty_hall_theorem)

        #Bridges of Königsberg
        bridges_point = Dot(axis.c2p(7,2.2))
        bridges_theorem = Text("Bridges of Königsberg").scale(0.5).next_to(bridges_point, DR, buff=0.05)
        bridges = VGroup(bridges_point, bridges_theorem)

        #Basel Problem
        basel_point = Dot(axis.c2p(3,4.5))
        basel_theorem = Text("Basel Problem").scale(0.5).next_to(basel_point, DR, buff=0.05)
        basel = VGroup(basel_point, basel_theorem)

        #Differentiation
        deriv_point = Dot(axis.c2p(10, 1))
        deriv_theorem = Text("Calc 1 Differentiation").scale(0.5).next_to(deriv_point, UR, buff=0.05)
        deriv = VGroup(deriv_point, deriv_theorem)

        #Integration
        integration_point = Dot(axis.c2p(5, 3))
        integration_theorem = Text("Calc 2 Integration").scale(0.5).next_to(integration_point, DR, buff=0.05)
        integration = VGroup(integration_point, integration_theorem)

        #Hairy Ball
        hairy_ball_point = Dot(axis.c2p(2, 6))
        hairy_ball_theorem = Text("Hairy Ball Theorem").scale(0.5).next_to(hairy_ball_point, UR, buff=0.05)
        hairy_ball = VGroup(hairy_ball_point, hairy_ball_theorem)

        problem_points = VGroup(
            fermat, arithmetic, ts, monty_hall, bridges, basel, deriv, integration, hairy_ball
        )
        self.play(
            LaggedStart(
                FadeIn(ts), FadeIn(monty_hall), FadeIn(bridges), FadeIn(basel),
                FadeIn(deriv), FadeIn(integration), FadeIn(hairy_ball)
            ),
            run_time=3
        )
        self.wait(3)

        self.play(axis.animate.shift(UP*1 + RIGHT*2))
        self.play(axis.animate.shift(DOWN*2 + LEFT*1))
        self.play(axis.animate.shift(UP*1 + LEFT*1))

        self.wait(2)
        self.play(
            hairy_ball.animate.shift(DOWN*0.5, RIGHT*0.2),
            integration.animate.shift(UP*1),
            basel.animate.shift(LEFT*1+UP*0.2),
            bridges.animate.shift(RIGHT*2+UP*1)
        )
        self.play(
            hairy_ball.animate.shift(UP*0.5, LEFT*0.2),
            integration.animate.shift(DOWN*1),
            basel.animate.shift(RIGHT*1+DOWN*0.2),
            bridges.animate.shift(LEFT*2+DOWN*1)
        )
        self.wait(4)
        self.play(fermat.animate.shift(UP*2))
        self.wait(2)
        self.play(arithmetic.animate.shift(UP*5.5+LEFT*0.5))
        self.wait(2)

        self.play(Create(dividers))
        self.wait()
        four_dots_point = Dot(axis.get_center()).scale(1.5).set_color(BLUE)
        
        self.play(FadeIn(four_dots_point))  
        self.play(FadeOut(problem_points))
        self.wait(3)
        four_dots_theorem = Text("Four Points, Two Distances").scale(0.6).next_to(four_dots_point, UL, buff=0.05).set_color(BLUE)
        self.play(FadeIn(four_dots_theorem))
        four_dots = VGroup(four_dots_theorem, four_dots_point)

        other_points = VGroup()

        for _ in range(7):
            other_points.add(
                Dot().move_to(four_dots_point.get_center() + np.array([rd.uniform(-0.4,0.4), rd.uniform(-0.4,0.4),0]))
            )
        self.play(FadeIn(other_points))
        

        self.wait()



class QuestionProposal(Scene):
    def construct(self):
        scale_factor = ValueTracker(2)
        square = Square(
            stroke_color = PURPLE_D,
            stroke_width = 5,
            side_length = scale_factor.get_value(),
        ).set(opacity=0).set_opacity(0)

        ur_dot = always_redraw(lambda: Dot(
            color = WHITE,
        ).move_to(square, UR).shift(UP*scale_factor.get_value()*0.03+RIGHT*scale_factor.get_value()*0.03)
        )
        ul_dot = always_redraw(lambda: Dot(
            color = WHITE,
        ).move_to(square, UL).shift(UP*scale_factor.get_value()*0.03+LEFT*scale_factor.get_value()*0.03)
        )
        dr_dot = always_redraw(lambda: Dot(
            color = WHITE,
        ).move_to(square, DR).shift(DOWN*scale_factor.get_value()*0.03+RIGHT*scale_factor.get_value()*0.03)
        )
        dl_dot = always_redraw(lambda: Dot(
            color = WHITE,
        ).move_to(square, DL).shift(DOWN*scale_factor.get_value()*0.03+LEFT*scale_factor.get_value()*0.03)
        )

        square_corners = VGroup(ur_dot, ul_dot, dr_dot, dl_dot)

        self.play(Create(square))

        self.wait(1.5)

        self.play(Create(square_corners))

        self.wait(2)

        square_corner_list = list(square_corners)
        sides = []
        for i in range(len(square_corner_list)):
            for j in range(i+1, len(square_corner_list)):
                sides.append(Group(square_corner_list[i], square_corner_list[j]))
        sides.pop(2)
        sides.pop(2)

        square_edges = list(draw_edges(square_corners, inner_lines=False))

        self.play(
            *[Create(mob) for mob in square_edges]
        )
                
        for side in sides:
            rect = SurroundingRectangle(side)
            self.play(Create(rect), run_time=0.5)
            self.play(Uncreate(rect.reverse_direction()), run_time=0.8)

        self.wait(2)

        ur_dl_line = Line(
            color = PURPLE_B,
            stroke_width = 5,
            start = ur_dot,
            end = dl_dot
        )

        ul_dr_line = Line(
            color = PURPLE_B,
            stroke_width = 5,
            start = ul_dot,
            end = dr_dot
        )

        self.play(Create(ur_dl_line), Create(ul_dr_line))
        
        diagonals = VGroup(ur_dl_line, ul_dr_line)

        square_group = VGroup(square, diagonals)
        line_group = VGroup(ur_dl_line, ul_dr_line)
        for edge in square_edges:
            line_group.add(edge)
        for edge in square_edges:
            square_group.add(edge)
        self.wait(4)
        self.play(square_group.animate.shift(LEFT*2))
        six_different_edges = Text("6 Different Edges").move_to([2,1,0]).scale(0.8)
        two_different_lengths = Text("2 Different Lengths").move_to([2.2,0.2,0]).scale(0.8)
        short_green = Text("- Short Outside Edges (Green)").move_to([3.5,-0.4,0]).scale(0.6).set_color(GREEN)
        long_red = Text("- Long Diagonals (Red)").move_to([3,-1,0]).scale(0.6).set_color(RED)
        edge_length_observations = VGroup(six_different_edges, two_different_lengths, short_green, long_red)
        self.play(Write(six_different_edges))
        self.wait()
        for mob in line_group:
            self.play(mob.animate.set_color(BLUE_B), run_time=0.5)
        self.wait(1.5)
        self.play(Write(two_different_lengths))
        self.wait()
        self.play(Write(short_green))
        for mob in square_edges:
            self.play(mob.animate.set_color(GREEN), run_time=0.5)
        self.play(Write(long_red))
        self.wait(1)
        for mob in diagonals:
            self.play(mob.animate.set_color(RED), run_time=0.5)
        self.wait(3)

        self.wait(2)
        self.play(
            edge_length_observations.animate.shift(DOWN*1),
            (mob.animate.scale(1.5).shift(DOWN*1+LEFT*1) for mob in [square, line_group]),
        )
        question = Text("""
                            How many, and in what ways can you arrange four dots such
                            that the connecting lines will have only two distances?""").scale(0.8).shift(UP*2)
        question_rect = SurroundingRectangle(question, color = PURPLE_E)
        question_group = VGroup(question, question_rect)
        self.wait()
        self.play(Write(question), run_time=3)
        self.play(Create(question_rect))
        self.wait(15)
        self.play(FadeOut(question_group))

        self.play(
            *[mob.animate.set_color(PURPLE_B) for mob in square_edges],
            *[mob.animate.set_color(PURPLE_B) for mob in diagonals]
        )
        
        self.play(
            FadeOut(edge_length_observations),
            square_group.animate.move_to([-3.5,1.5,0]).set_width(1.6)
        )

        self.wait(2)

        #THE OTHER 5 WAYS
        rhombus_unorganized = VGroup(
            Dot().shift(RIGHT*0.3),
            Dot().shift(LEFT*1.2+UP*0.4),
            Dot().shift(DOWN*0.3+LEFT*0.1),
            Dot().shift(UP*0.7+RIGHT*0.8)
        ).move_to([0,1.5,0])

        kite_unorganized = VGroup(
            Dot().shift(RIGHT*0.5),
            Dot().shift(LEFT*0.6+UP*0.9),
            Dot().shift(DOWN*1.2+LEFT*0.3),
            Dot().shift(UP*0.3+RIGHT*0.5)
        ).move_to([3.5,1.5,0])

        isosceles_unorganized = VGroup(
            Dot().shift(RIGHT*0.5),
            Dot().shift(LEFT*0.7+UP*0.3),
            Dot().shift(DOWN*0.9+LEFT*0.3),
            Dot().shift(UP*0.1+RIGHT*0.8)
        ).move_to([-3.5,-1.5,0])

        centered_triangle_unorganized = VGroup(
            Dot().shift(RIGHT*0.5),
            Dot().shift(LEFT*0.6+UP*0.4),
            Dot().shift(DOWN*1.1+LEFT*0.4),
            Dot().shift(UP*0.1+RIGHT*0.3)
        ).move_to([0,-1.5,0])

        trapezoid_unorganized = VGroup(
            Dot().shift(RIGHT*0.9),
            Dot().shift(LEFT*0.6+UP*0.3),
            Dot().shift(DOWN*0.4+LEFT*0.3),
            Dot().shift(UP*0.4+RIGHT*0.1)
        ).move_to([3.5,-1.5,0])

        unorganized_dots = VGroup(
            rhombus_unorganized,
            kite_unorganized,
            isosceles_unorganized,
            centered_triangle_unorganized,
            trapezoid_unorganized
        )
        self.play(
            FadeIn(unorganized_dots)
        )
        self.wait(3)

        #KITE
        kite_edges = draw_edges(kite_unorganized, update_edges=True)
        self.play(FadeIn(kite_edges))
        self.play(
            *[
                mob.animate.move_to(2*kite_array[i] + kite_unorganized.get_center()) for i, mob in enumerate(kite_unorganized)
            ]
        )
        rhombus_edges = draw_edges(rhombus_unorganized, update_edges=True)
        self.play(FadeIn(rhombus_edges))
        self.play(
            *[
                mob.animate.move_to(rhombus_array[i] + rhombus_unorganized.get_center()) for i, mob in enumerate(rhombus_unorganized)
            ]
        )
        centered_triangle_edges = draw_edges(centered_triangle_unorganized, update_edges=True)
        self.play(FadeIn(centered_triangle_edges))
        self.play(
            *[
                mob.animate.move_to(2*centered_triangle_array[i] + centered_triangle_unorganized.get_center()) for i, mob in enumerate(centered_triangle_unorganized)
            ]
        )
        isosceles_edges = draw_edges(isosceles_unorganized, update_edges=True)
        self.play(FadeIn(isosceles_edges))
        self.play(
            *[
                mob.animate.move_to(2*isosceles_array[i] + isosceles_unorganized.get_center()+UP*0.8) for i, mob in enumerate(isosceles_unorganized)
            ]
        )

        self.wait(4)

        self.play(
            Circumscribe(trapezoid_unorganized, Circle)
        )
        self.wait(3)
        trapezoid_edges = draw_edges(trapezoid_unorganized, update_edges=True)
        self.play(FadeIn(trapezoid_edges))
        self.play(
            *[
                mob.animate.move_to(2*trapezoid_array[i] + trapezoid_unorganized.get_center()) for i, mob in enumerate(trapezoid_unorganized)
            ]
        )
        self.wait(2)
        """"
        Numbers = MathTex("0", "1", "2", "3", "4", "5")

        self.play(
            *[
                FadeIn(mob.move_to(trapezoid_edges[i].get_center()+UP*0.2)) for i, mob in enumerate(Numbers)
            ],
        )
        self.play(
            *[
                FadeIn(mob.move_to(rhombus_edges[i].get_center()+UP*0.2)) for i, mob in enumerate(Numbers)
            ],
        )
        self.play(
            *[
                FadeIn(mob.move_to(kite_edges[i].get_center()+UP*0.2)) for i, mob in enumerate(Numbers)
            ],
        )
        self.play(
            *[
                FadeIn(mob.move_to(centered_triangle_edges[i].get_center()+UP*0.2)) for i, mob in enumerate(Numbers)
            ],
        )
        self.play(
            *[
                FadeIn(mob.move_to(isosceles_edges[i].get_center()+UP*0.2)) for i, mob in enumerate(Numbers)
            ],
        )
        """

        short_edges = Group(
            trapezoid_edges[0], trapezoid_edges[1], trapezoid_edges[4],
            rhombus_edges[0], rhombus_edges[1], rhombus_edges[3], rhombus_edges[4], rhombus_edges[5],
            kite_edges[1], kite_edges[2],
            centered_triangle_edges[0], centered_triangle_edges[1], centered_triangle_edges[2],
            isosceles_edges[0], isosceles_edges[3], isosceles_edges[4], isosceles_edges[5], 
            square_edges[0], square_edges[1], square_edges[2], square_edges[3]
        )
        long_edges = Group(
            trapezoid_edges[2], trapezoid_edges[3], trapezoid_edges[5],
            rhombus_edges[2],
            kite_edges[0], kite_edges[3], kite_edges[4], kite_edges[5],
            centered_triangle_edges[3], centered_triangle_edges[4], centered_triangle_edges[5],
            isosceles_edges[1], isosceles_edges[2],
            diagonals[0], diagonals[1]
        )

        self.play(
            *[
                mob.animate.set_color(GREEN) for mob in short_edges
            ]
        )
        self.play(
            *[
                mob.animate.set_color(RED) for mob in long_edges
            ]
        )
        self.wait(3)

        #SCENE TRANSITIONS
"""
        IDEA I HAD ORIGINALLY, FOUND A BETTER PROOF IN THE PROCESS OF ANIMATING
        
        long_symbol = MathTex("L\\text{: Long Edge}").scale(2)
        short_symbol = MathTex("S\\text{: Short Edge}").scale(2)

        self.wait()

        self.play(FadeIn(long_symbol))
        self.play(
            FadeIn(short_symbol.shift(DOWN*1)),
            long_symbol.animate.shift(UP*1)
        )
        self.wait(3)

        fade_all(self)

        possibilities = Group(
            MathTex("AAAAAB"),
            MathTex("AAAABB"),
            MathTex("AAABBB"),
        )
        self.play(*[Write(combination.shift((-2+i)*DOWN)) for i, combination in enumerate(possibilities)])
        self.wait(2)

        temp = possibilities[0].copy()
        self.add(temp)
        self.remove(possibilities[0])
        for i in range(1, len(possibilities)):
            self.play(
                FadeOut(temp),
                *[FadeIn(combination) for combination in possibilities]
            )

            temp = possibilities[i].copy()
            self.play(
                FadeOut(possibilities),
                temp.animate.move_to([0,2,0])
            )
            self.wait(2)

            match i:
                case 1:
                    None
            self.wait()
        self.play(FadeOut(possibilities[1:6]))
        self.wait(2)

"""
class Methodology(Scene):
    def construct(self):

        rhombus_origin = Dot().set_opacity(0)
        rhombus_height = ValueTracker(2)

        Lemma1 = Text(
            """
            Three distinct points fully determine a triangle
            """
        ).move_to([0,2.5,0])
        self.play(Create(Lemma1), run_time=3)
        self.wait()

        #TRIANGLE ANIMATION
        triangle_corners = VGroup(
            Dot().shift(UP*1+RIGHT*1.2),
            Dot().shift(LEFT*1.2 + DOWN*0.3),
            Dot().shift(UP*0.5+LEFT*0.8),
        )
        triangle_edges = draw_edges(triangle_corners, update_edges=True)

        self.play(FadeIn(triangle_corners))
        self.play(FadeIn(triangle_edges))
        self.wait()
        self.play(triangle_corners[0].animate.shift(UP*0.7+RIGHT*0.4))
        self.wait()
        #TEMPORARY 4th DOT DEMONSTRATION - ALL QUADRILATERALS CAN BE FORMED FROM A TRIANGLE
        self.wait(2)
        fourth_dot = Dot().move_to([1,-2,0])
        
        
        self.play(
            FadeIn(fourth_dot)
        )
        self.wait(2)

        triangle_corners.add(fourth_dot)
        temp_edges = draw_edges(triangle_corners, update_edges=False)
        self.play(
            FadeIn(temp_edges)
        )
        self.wait(10)
        self.play(
            FadeOut(temp_edges),
            FadeOut(fourth_dot)
        )

        triangle_corners.remove(fourth_dot)
        triangle_edges = draw_edges(triangle_corners, update_edges=True)

        #EQUILATERAL CASE
        triangle_side_length = 3
        self.wait(3)
        self.play(
            *[
                mob.animate.move_to(triangle_side_length*equilateral_array[i]+[-0.5,-0.5,0]) for i, mob in enumerate(triangle_corners)
            ]
        )

        self.wait(2)
        self.play(Rotate(triangle_corners,3*PI/2))
        self.wait(2)
        self.play(triangle_corners[2].animate.shift(LEFT*1.5))
        self.wait(6)
        self.play(triangle_corners[2].animate.shift(RIGHT*1.5))
        self.wait(8)

        fourth_dot = Dot(color = BLUE).move_to((triangle_corners[0].get_center()+triangle_corners[1].get_center())/2 + RIGHT*1)
        
        bisecting_line = DashedLine(
            start = fourth_dot.get_center() + RIGHT*3,
            end = triangle_corners[2].get_center() + LEFT*4
        ).set_opacity(0.3)
        self.play(FadeIn(bisecting_line))
        self.wait(5)
        self.play(FadeIn(fourth_dot))
        triangle_corners.add(fourth_dot)
        
        #JUSTIFICATION FOR POINT LYING ON BISECTING LINE
        new_edges = Group()
        for i in range(len(triangle_corners)):
            new_edges.add(
                draw_edges(
                    Group(fourth_dot, triangle_corners[i]),
                    color = WHITE,
                    update_edges = True
                )
            )
        self.wait(3)
        self.play(FadeIn(new_edges))
        self.wait(3)

        self.play(
            (mob.animate.set_color(PURE_RED) for mob in new_edges[0:2])
        )
        self.wait(3)
        self.play(
            new_edges[2].animate.set_color(BLUE)
        )
        self.wait(3)
        self.play(fourth_dot.animate.shift(UP*0.5+RIGHT*0.5), run_time=2)
        self.wait(2)
        self.play(new_edges[1].animate.set_color(ORANGE), run_time=0.5)
        self.wait()
        self.play(new_edges[1].animate.set_color(PURE_RED))
        self.play(fourth_dot.animate.shift(DOWN*0.5+LEFT*0.5), run_time=2)
        self.wait(2)
        self.play(FadeOut(Lemma1))
        
        #PROOF OF 4/6 QUADRILATERAL
        self.wait()
        l = MathTex("l").next_to(triangle_edges[1].get_end(), UP, buff=-0.3).set_color(PURPLE)
        a = always_redraw(lambda: MathTex("a").next_to(new_edges[0].get_center(), UP, buff=0.2).set_color(PURE_RED))
        b = always_redraw(lambda: MathTex("b").next_to(new_edges[2].get_center(), UP, buff=0.15).set_color(BLUE))
        side_lengths = VGroup(l, a, b)
        self.play(FadeIn(side_lengths))
        self.wait(3)

        
        l_is_constant = MathTex("l\\text{ is (effectively) constant:}").scale(1).move_to([-3.5,2.5,0]).set_color(PURPLE)
        self.play(Create(l_is_constant))
        self.wait(5)
        first_case = MathTex("\\text{Case 1: }","a", "=", "l").move_to([4,2.4,0]).set_color_by_tex('l', PURPLE).set_color_by_tex('a', RED)
        first_case[0].set_color(WHITE)
        second_case = MathTex("\\text{Case 2: }", "b", "=", "l").move_to([4,1.7,0]).set_color_by_tex('b', BLUE).set_color_by_tex('l', PURPLE)
        third_case = MathTex("\\text{Case 3: }", "a", "=", "b").move_to([4,1,0]).set_color_by_tex('a', RED).set_color_by_tex('b', BLUE)
        third_case[0].set_color(WHITE)
        self.play(FadeIn(first_case))
        self.wait()
        self.play(FadeIn(second_case))
        self.wait()
        self.play(FadeIn(third_case))
        self.wait(6)
        a_lower = always_redraw(lambda: MathTex("a").next_to(new_edges[1].get_center(), DOWN, buff=0.2).set_color(PURE_RED))
        self.play(FadeIn(a_lower))

        self.play(
            fourth_dot.animate.shift(RIGHT*2.5),
            run_time = 4
        )
        self.wait()
        self.play(fourth_dot.animate.shift(RIGHT*1), run_time=4)
        self.wait(4)

        first_eq_pt = triangle_corners[2].get_center() + RIGHT*triangle_side_length*np.sqrt(3)
        second_eq_pt = triangle_corners[2].get_center() + RIGHT*triangle_side_length
        third_eq_pt = triangle_corners[2].get_center() + RIGHT*triangle_side_length*np.sqrt(3)/3
        fourth_eq_pt = triangle_corners[2].get_center() + LEFT*triangle_side_length

        #RHOMBUS
        self.play(
            fourth_dot.animate.move_to(first_eq_pt), run_time=3
        )
        self.play(fourth_dot.animate.set_color(PURE_RED))
        self.play(
            Group(
                triangle_edges,
                new_edges[0:2]
            ).animate.set_color(PURE_GREEN)
        )
        self.play(
            Circumscribe(first_case)
        )
        self.play(
            triangle_edges.animate.set_color(PURPLE),
            new_edges[0:2].animate.set_color(PURE_RED),
        )
        self.wait()
        self.play(fourth_dot.animate.set_color(WHITE))

        #KITE
        self.play(
            fourth_dot.animate.move_to(second_eq_pt), run_time=3
        )
        self.play(fourth_dot.animate.set_color(PURE_RED))
        self.play(
            Group(
                triangle_edges,
                new_edges[2]
            ).animate.set_color(PURE_GREEN)
        )
        self.play(
            Circumscribe(second_case)
        )
        self.play(
            triangle_edges.animate.set_color(PURPLE),
            new_edges[2].animate.set_color(BLUE),
        )
        self.wait()
        self.play(fourth_dot.animate.set_color(WHITE))

        #CENTERED EQUILATERAL
        self.play(
            fourth_dot.animate.move_to(third_eq_pt), run_time=3
        )
        self.play(fourth_dot.animate.set_color(PURE_RED))
        self.play(
            Group(
                new_edges
            ).animate.set_color(PURE_GREEN)
        )
        self.play(
            Circumscribe(third_case)
        )
        self.play(
            new_edges[0:2].animate.set_color(PURE_RED),
            new_edges[2].animate.set_color(BLUE),
        )
        self.wait()
        self.play(fourth_dot.animate.set_color(WHITE))
        #EQUILATERAL + ISOSCELES
        self.play(
            fourth_dot.animate.move_to(fourth_eq_pt), run_time=3
        )
        self.play(fourth_dot.animate.set_color(PURE_RED))
        self.play(
            Group(
                triangle_edges,
                new_edges[2]
            ).animate.set_color(PURE_GREEN)
        )
        self.play(
            Circumscribe(second_case)
        )
        self.play(
            triangle_edges.animate.set_color(PURPLE),
            new_edges[2].animate.set_color(BLUE),
        )
        self.wait()
        self.play(fourth_dot.animate.set_color(WHITE))
        self.wait(3)
        self.play(fourth_dot.animate.shift(LEFT*1), run_time=4)

        self.wait(2)


class PentagonProof(Scene):
    def construct(self):
        pentagon_corners = VGroup()
        for i in range(5):
            pentagon_corners.add(
                Dot().move_to(2*pentagon_array[i])
            )
        
        self.play(FadeIn(pentagon_corners))
        self.wait(3)
        pentagon_edges = draw_edges(pentagon_corners, update_edges=True)
        self.play(Write(pentagon_edges), run_time=2)
        self.wait(2)
        self.play(Rotate(pentagon_corners, angle=PI/5))
        self.wait(2)
        trapezoid_edges = VGroup(pentagon_edges[0], pentagon_edges[2:4], pentagon_edges[5:7], pentagon_edges[9])
        self.play(*[FadeOut(mob) for mob in [pentagon_edges[1], pentagon_edges[4], pentagon_edges[7:9], pentagon_corners[2]]])
        self.wait(2)
        self.play(mob.animate.set_color(PURE_BLUE) for mob in (pentagon_edges[2], pentagon_edges[5:7]))
        self.wait(2)
        self.play(mob.animate.set_color(PURE_GREEN) for mob in (pentagon_edges[0], pentagon_edges[3], pentagon_edges[9]))
        self.wait(2)
        self.play(trapezoid_edges.animate.set_color(PURPLE))
        self.wait()

        







