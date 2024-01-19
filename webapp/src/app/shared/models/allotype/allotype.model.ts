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
export class HiResAllele {
    typing : string;
    ciwd : string;
    exon3Motif : string;
    utr3Motif : string;
    experimental : boolean;
    exprLevel : string;
    tce : string;
}

export class Allotype {

    typing : string;
    submittedHlaB: string;
    exprLevel : string;
    exprAnnotType : string;
    tce : string;
    alleles : HiResAllele[];

    initiatedCall: boolean;
    sharedIndex: number;
    highlighted: boolean;
    
    constructor(typing: string) {
        this.typing = typing;
        this.initiatedCall = false;
        this.highlighted = false;
    }

}