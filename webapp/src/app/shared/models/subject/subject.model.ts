/*
 * Copyright (c) 2024 NMDP.
 *
 * This file is part of the NMDP DP Tool 
 * (see https://github.com/nmdp-bioinformatics/dpb1-expression).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
import { Allotype } from '@app/shared/models/allotype/allotype.model'

export class Subject {

    allotypes : Allotype[];
    type: string;
    id: string;
    index: number;

    genotype : string = null;
    genotypeDisplayed : string[] = ['', ''];
    
    expressionMatch: string[] = ["","",""];
    hlaMatch: string[];
    sharedAllotype : string = null;
    // flippedPatient: boolean = false;
    flipped: boolean = false;
    loading: boolean = false;
    rank: number = null;
    annotated : boolean = false;
    error : string;

    directionality : string;
    matchGrade : string;
    tceMatch : string;
    mismatchedAllelePatExprLevel : string;

    constructor(hlaBallotypes: string[], type: string, id : string, index: number) {
        this.allotypes = hlaBallotypes.map(allo => new Allotype(allo));
        this.type = type;
        this.id = id;
        this.index = index;
    }

    setHlaBAllotypes(hlaBallotypes: string[]){
        for (var i = 0; i < this.allotypes.length; i++){
            this.allotypes[i].typing = hlaBallotypes[i];
            this.genotypeDisplayed[i] = 'DPB1*' + hlaBallotypes[i];
        }
    }

}