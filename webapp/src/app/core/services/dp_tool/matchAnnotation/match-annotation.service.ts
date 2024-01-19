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
import { Injectable } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model';
import { HttpClient } from '@angular/common/http';
import { environment } from "../../../../../environments/environment";
import { BehaviorSubject } from 'rxjs';

@Injectable({
  providedIn: 'root'
})
export class MatchAnnotationService {
  private baseURL = environment.apiUrl + '/sorting/ids';
  private _matchingSource = new BehaviorSubject<void>(null);
  matching$ = this._matchingSource.asObservable()
  
  // private _flippingSource = new BehaviorSubject<string>(null);
  // flipping$ = this._flippingSource.asObservable();

  constructor(private httpClient: HttpClient) { }

  changeMatching() {
    this._matchingSource.next();
  }

  // changeFlipping(id : string) {
  //   this._flippingSource.next(id);
  // }

  getMatchInfo(patient: Subject, donors: Subject[], methods : string[] = ['expression', 'tce']): Promise<any> {
    donors = donors.filter(d => !d.error);
    const params = this._formatInput(patient, donors, methods);
    // console.log("loading");
    donors.forEach(donor => donor.loading = true)
    return this.httpClient.post(this.baseURL, params)
      .toPromise()
      .then((response: Object[]) => {
        response.forEach(res => {
          if (res.hasOwnProperty('message')){
            // const id = res['id'];
            // const donorsFiltered = donors.filter(d => d.id == id);
            // if (donorsFiltered.length != 1){
            //   console.log("Check if IDs are made unique and properly submitted and received.")
            // }
            // const donor = donorsFiltered[0];
            // donor.error = res['message'];
          } else {
            const genoRecip = res['genotype_patient']['genotype'].split('+')
            let sharedIndex = null;
            if (res['grade'] == 'AA'){
              sharedIndex = 2;
            } else {
              sharedIndex = genoRecip.indexOf(res['matched_alleles_pat'][0]);
            }
            const id = res['genotype_donor']['id'];
            const donorsFiltered = donors.filter(d => d.id == id);
            if (donorsFiltered.length != 1){
              console.log("Check if IDs are made unique and properly submitted and received.")
            }
            const donor = donorsFiltered[0];
            donor.mismatchedAllelePatExprLevel = res['annotation']['expr_match']
            donor.directionality = res['directionality'];
            donor.hlaMatch = res['grade'].split('');
            donor.tceMatch = res['annotation']['tce_match'];
            patient.genotype = res['genotype_patient']['genotype'];
            donor.genotype = res['genotype_donor']['genotype'];
            donor.annotated = true;
            donor.loading = false;
            if (res.hasOwnProperty('rank') && donors.length > 1){
              donor.rank = res['rank'];
            } else {
              donor.rank = null;
            }
            donor.error = null;
            // if ((donor.genotype.indexOf('DPB1*14:01') > 0) &&
            //     (res['genotype_recip']['genotype'] == 'DPB1*14:01+DPB1*03:01')){
            // }
            // if (res['genotype_donor_flipped']){
            //   this.changeFlipping(id);
            // }
            donor.flipped = res['genotype_donor_flipped'];
            this._assignAlleleInfo(donor, 'donor', res, sharedIndex);
            this._assignAlleleInfo(patient, 'patient', res);
          }
        })
      })
      .catch((err) => {
        // if (donors.length == 1){
        donors.forEach(donor => {
          // let donor = donors[0]
          donor.loading = false;
          donor.error = err['error'];
          donor.rank = null;
          console.log('errored!', donor);
          console.log(err);
          console.log(params);
        })
        // }
      })
  }

  private _assignAlleleInfo(subject : Subject, type : string, res : Object,
      sharedIndex : number = null){
    subject.allotypes.forEach((allo, index) => {
      if (subject.flipped){
        index = 1 - index;
      }
      const dpb1 = 'DPB1*' + allo.typing
      if (sharedIndex != null && res['matched_alleles_don'].indexOf(dpb1) >= 0){
        allo.sharedIndex = sharedIndex;
      } else {
        allo.sharedIndex = null;
      }
      let allele_index = 'allotype_' + ((0 == index) ? 'one' : 'two');
      let allele = res['genotype_' + type][allele_index];
      allo.submittedHlaB = allele['typing'].replace('DPB1*', '');
      allo.exprLevel = allele['feats']['anns']['expression'];
      allo.exprAnnotType = allele['feats']['anns']['expression_experimental'];
      allo.tce = allele['feats']['anns']['tce'];
      // let header = '';
      // if (allele['alleles_hi_res']){
      //   header = 'alleles_hi_res'
      // } else if (allele['alleles'] && (allo.exprLevel.includes('?') || allo.exprLevel.includes('~'))) {
      //   header = 'alleles'
      // }
      console.log(allele);
      if (allele['feats_expanded']) {
        allo.alleles = allele['feats_expanded'].map(a => {
          return {'typing' : a['allele_two_field'],
                  'ciwd' : a['CIWD_TOTAL_combined'],
                  'experimental' : a['expression_experimental'],
                  'exprLevel' : a['expression'],
                  'tce' : a['tce']}
        });
      }
    })
  }

  private _formatInput(patient: Subject, donors: Subject[], methods : string[]) {
    let input = {"dpb1_genotype_patient" : this._formatGenotypes([patient])[0]['genotype'],
            "dpb1_genotype_donors"  : this._formatGenotypes(donors),
            "sorting_methods" : methods};
    return input
  }

  private _formatGenotypes(subjects: Subject[]){
    return subjects.filter(sub => sub.genotypeDisplayed.join('') != '')
        .map(sub => {
                let genotype = sub.genotypeDisplayed.join('+') //.map(allo => 'DPB1*' + allo)
                return {'id' : sub.id, 'genotype' : genotype}
    });
  }

}
