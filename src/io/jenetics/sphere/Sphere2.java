package io.jenetics.sphere;

import java.util.List;
import java.util.function.Function;
import static java.util.stream.Collectors.toList;


import io.jenetics.Chromosome;
import io.jenetics.Genotype;
import io.jenetics.Mutator;
import io.jenetics.engine.Codec;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.util.ISeq;
import io.jenetics.util.RandomRegistry;


import io.jenetics.ext.SingleNodeCrossover;


import io.jenetics.prog.ProgramChromosome;
import io.jenetics.prog.ProgramGene;
import io.jenetics.prog.op.EphemeralConst;
import io.jenetics.prog.op.MathOp;
import io.jenetics.prog.op.Op;
import io.jenetics.prog.op.Var;
/**
 * @ClassName Sphere2
 * @Description: TODO
 * @Author: lichenxingbeijing@163.com
 */
public class Sphere2 {
	public static final double XMIN = -10.0;
	public static final double XMAX =  10.0;
    public static final Op<Double> PDIV = Op.of("pdiv", 2, v -> (v[1] == 0) ? 1 : v[0]/v[1]);


    // Definition of the allowed operations.
    static final ISeq<Op<Double>> OPERATIONS = ISeq.of(
            MathOp.ADD,
            MathOp.SUB,
            MathOp.MUL,
            PDIV
    );


    // Definition of the terminals.
    static final ISeq<Op<Double>> TERMINALS = ISeq.of(
            EphemeralConst.of(() -> RandomRegistry.random().nextDouble())
    );


    static final Codec<ISeq<Function<Double[], Double>>, ProgramGene<Double>> CODEC =
            Codec.of(
                    Genotype.of(
                            // First 'program'
                            ProgramChromosome.of(
                                    5,
                                    ch -> ch.root().size() <= 50,
                                    OPERATIONS,
                                    TERMINALS
                            ),
                            // Second 'program'
                            ProgramChromosome.of(
                                    5,
                                    ch -> ch.root().size() <= 50,
                                    OPERATIONS,
                                    TERMINALS
                            )
                    ),
                    gt -> gt.stream()
                            .map(c -> (ProgramGene<Double>) c.gene())
                            .collect(ISeq.toISeq())
            );


    static double fitness(final ISeq<Function<Double[], Double>> programs) {
        assert programs.size() == 2;
        return programs.stream()
                .map(p -> p.apply(new Double[0]))
                .map(x -> Math.max(x, XMIN))
                .map(x -> Math.min(x, XMAX))
                .reduce(0.0, (acc, x) -> acc + (x * x));
    }


    public static void main(final String[] args) {
    	long startTime = System.currentTimeMillis();
        final Engine<ProgramGene<Double>, Double> engine = Engine
                .builder(Sphere2::fitness, CODEC)
                .minimizing()
                .alterers(
                        new SingleNodeCrossover<>(),
                        new Mutator<>())
                .build();


        final Genotype<ProgramGene<Double>> gt = engine.stream()
                .limit(4000)
                .collect(EvolutionResult.toBestGenotype());


        final ISeq<Function<Double[], Double>> programs = CODEC.decode(gt);

        List<Double> best = programs.stream()
        		.map(p -> p.apply(new Double[0]))
        		.map(x -> Math.max(x, XMIN))
        		.map(x -> Math.min(x,XMAX))
        		.collect(toList());

        System.out.format("%.6f\n", fitness(programs));
        System.out.println(best.stream().map(x -> String.format("%.6f", x)).collect(toList()));
        long endTime = System.currentTimeMillis();
        System.out.println(endTime-startTime + " ms");
    }
}


